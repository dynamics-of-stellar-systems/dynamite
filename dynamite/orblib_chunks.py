"""Merging orbit libraries that were integrated in chunks.

When an orbit library is built by several processes covering disjoint orbit
ranges, each process writes its own output files. This module stitches those
back into the single files the rest of DYNAMITE expects.

Both orbit library streams are sequences of Fortran unformatted records: a
fixed header describing the library, then a body with one group of records per
orbit (``_qgrid.dat``) or per orbit and aperture (``_losvd_hist.dat``), then a
one-byte record closing the file. Chunks cover consecutive, disjoint orbit
ranges written in the same format, so merging is a byte-level operation: keep
the first chunk's header, concatenate every chunk's body in order, and append a
single closing record.

Merged files are byte-identical to the output of a single-process run.
"""

import os
import shutil
import struct

import numpy as np

# number of Fortran records making up each file's header
HEADER_RECORDS = {'qgrid': 5, 'losvd_hist': 1}


def header_length(buf, n_records):
    """Byte offset just past the first ``n_records`` Fortran records.

    Parameters
    ----------
    buf : 1d numpy array of uint8
        the file contents
    n_records : int
        number of records making up the header

    Returns
    -------
    int
        offset of the first byte after the header

    Raises
    ------
    ValueError
        if the file is truncated, or a record's leading and trailing length
        markers disagree, which means it is not the expected format.

    """
    pos = 0
    for _ in range(n_records):
        if pos + 4 > buf.nbytes:
            raise ValueError('truncated file while reading the header')
        n = int(np.frombuffer(buf[pos:pos + 4], dtype=np.int32)[0])
        trail = int(np.frombuffer(buf[pos + 4 + n:pos + 8 + n],
                                  dtype=np.int32)[0])
        if trail != n:
            raise ValueError('record length markers disagree '
                             f'({n} vs {trail}) - not an orbit library file?')
        pos += 8 + n
    return pos


def footer_offset(buf):
    """Byte offset of the one-byte record that closes an orbit library file.

    Every chunk writes this record, so a merged file must carry exactly one of
    them rather than one per chunk. It is located from the file's trailing
    length marker, so this costs nothing on a large file.

    Parameters
    ----------
    buf : 1d numpy array of uint8
        the file contents

    Returns
    -------
    int
        offset of the first byte of the closing record

    Raises
    ------
    ValueError
        if the closing record's length markers disagree.

    """
    n = int(np.frombuffer(buf[-4:], dtype=np.int32)[0])
    start = buf.nbytes - 8 - n
    lead = int(np.frombuffer(buf[start:start + 4], dtype=np.int32)[0])
    if lead != n:
        raise ValueError('trailing record length markers disagree '
                         f'({lead} vs {n}) - file may be incomplete')
    return start


def merge_files(chunk_files, out_file, kind, n_orbits=None):
    """Concatenate chunk files of one kind into a single orbit library file.

    Parameters
    ----------
    chunk_files : list of string
        the chunks' file names, **in ascending orbit order**
    out_file : string
        file to write
    kind : string
        ``'qgrid'`` or ``'losvd_hist'``, selecting the header length
    n_orbits : int, optional
        total number of orbits. Written into the ``qgrid`` header, which
        records it; ignored for other kinds. The default is None.

    Returns
    -------
    string
        ``out_file``

    """
    bufs = [np.fromfile(f, dtype=np.uint8) for f in chunk_files]
    head_len = header_length(bufs[0], HEADER_RECORDS[kind])
    header = bytearray(bufs[0][:head_len].tobytes())
    if kind == 'qgrid' and n_orbits is not None:
        # the first record's payload begins at byte 4 and starts with n_orbits
        struct.pack_into('<i', header, 4, n_orbits)
    with open(out_file, 'wb') as f:
        f.write(header)
        for buf in bufs:
            f.write(buf[header_length(buf, HEADER_RECORDS[kind]):
                        footer_offset(buf)].tobytes())
        f.write(bufs[-1][footer_offset(bufs[-1]):].tobytes())
    return out_file


def merge_orbclass(datfil, fileroot, chunk_tags):
    """Concatenate one orbit family's ``orbclass.out`` chunks.

    Unlike the two binary streams this one is text and has no header: the line
    ``integrator_setup_write`` writes to unit 30 is discarded, because
    ``output_setup`` reopens unit 30 with ``status="replace"`` afterwards. So
    the file is nothing but one ``integrator_moments`` block per orbit, and
    merging is a concatenation in ascending orbit order.

    Parameters
    ----------
    datfil : string
        the model's ``datfil`` directory
    fileroot : string
        ``'orblib'`` or ``'orblibbox'``
    chunk_tags : list of string
        per-chunk file name suffixes, in ascending orbit order

    Returns
    -------
    string
        the merged file name

    Raises
    ------
    FileNotFoundError
        if any chunk file is missing, naming the first one.

    """
    parts = [os.path.join(datfil, f'{fileroot}{tag}.dat_orbclass.out')
             for tag in chunk_tags]
    for p in parts:
        if not os.path.isfile(p):
            raise FileNotFoundError(f'missing orbit class chunk: {p}')
    out_file = os.path.join(datfil, f'{fileroot}.dat_orbclass.out')
    # write to a temporary name and rename, so a concurrently evaluated model
    # sharing this library never sees a half-written file
    tmp = f'{out_file}.tmp{os.getpid()}'
    with open(tmp, 'wb') as out:
        for p in parts:
            with open(p, 'rb') as f:
                shutil.copyfileobj(f, out)
    os.replace(tmp, out_file)
    return out_file


def remove_chunks(datfil, fileroot, chunk_tags,
                  kinds=('qgrid', 'losvd_hist')):
    """Delete one orbit family's chunk files, ignoring any already gone.

    Used to discard a process' chunks when another process has already merged
    the library, so that concurrently evaluated models sharing an orbit library
    do not leave a full extra copy of it on disk.

    Parameters
    ----------
    datfil : string
        the model's ``datfil`` directory
    fileroot : string
        ``'orblib'`` or ``'orblibbox'``
    chunk_tags : list of string
        per-chunk file name suffixes
    kinds : tuple of string, optional
        which output streams to delete. The default is
        ``('qgrid', 'losvd_hist')``.

    """
    for kind in kinds:
        for tag in chunk_tags:
            for suffix in ('', '.tmp'):
                path = os.path.join(
                    datfil, f'{fileroot}{tag}_{kind}.dat{suffix}')
                if os.path.isfile(path):
                    os.remove(path)
    for tag in chunk_tags:
        path = os.path.join(datfil, f'{fileroot}{tag}.dat_orbclass.out')
        if os.path.isfile(path):
            os.remove(path)


def merge_chunks(datfil, fileroot, chunk_tags, n_orbits,
                 kinds=('qgrid', 'losvd_hist'), cleanup=True, out_tag=''):
    """Merge one orbit family's chunk files, and delete the chunks.

    Parameters
    ----------
    datfil : string
        the model's ``datfil`` directory
    fileroot : string
        ``'orblib'`` or ``'orblibbox'``
    chunk_tags : list of string
        per-chunk file name suffixes, in ascending orbit order, as used when
        the chunks' input files were written
    n_orbits : int
        total number of orbits in the merged library
    kinds : tuple of string, optional
        which output streams to merge. The default is
        ``('qgrid', 'losvd_hist')``.
    cleanup : bool, optional
        whether to delete the chunk files once merged. The default is True.
    out_tag : string, optional
        suffix for the merged file names. Concurrently evaluated models sharing
        an orbit library would otherwise merge into the same file at the same
        time; giving each its own name and renaming afterwards keeps that safe.
        The default is '', i.e. the plain names.

    Returns
    -------
    list of string
        the merged file names

    Raises
    ------
    FileNotFoundError
        if any chunk file is missing, naming the first one.

    """
    merged = []
    for kind in kinds:
        parts = [os.path.join(datfil, f'{fileroot}{tag}_{kind}.dat')
                 for tag in chunk_tags]
        for p in parts:
            if not os.path.isfile(p):
                raise FileNotFoundError(f'missing orbit library chunk: {p}')
        merged.append(merge_files(
            parts, os.path.join(datfil, f'{fileroot}_{kind}.dat{out_tag}'),
            kind, n_orbits))
    # not returned with the others: this one is text, is not compressed, and is
    # already at its final name, so the caller must not post-process it
    merge_orbclass(datfil, fileroot, chunk_tags)
    if cleanup:
        remove_chunks(datfil, fileroot, chunk_tags, kinds)
    return merged
