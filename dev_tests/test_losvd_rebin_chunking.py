"""Check the chunked BLAS rebin matches the einsum it replaced.

Mirrors the contraction in Histogram.rebin_orblib_to_observations, including
a chunk size small enough to exercise a ragged final chunk. Run this file.
"""
import numpy as np


def reference(y, f):
    return np.einsum('ijk,lj->ilk', y, f, optimize=False)


def chunked(y, f, chunk):
    n_orb, n_vbin, n_ap = y.shape
    out = np.empty((n_orb, f.shape[0], n_ap), dtype=y.dtype)
    for start in range(0, n_orb, chunk):
        end = min(start + chunk, n_orb)
        blk = y[start:end].transpose(1, 0, 2).reshape(n_vbin, -1)
        out[start:end] = (f @ blk).reshape(f.shape[0], end - start,
                                           n_ap).transpose(1, 0, 2)
    return out


def demo():
    rng = np.random.default_rng(0)
    for n_orb, n_vbin, n_ap, n_data in ((37, 20, 9, 13), (5, 4, 1, 2)):
        y = rng.random((n_orb, n_vbin, n_ap))
        f = rng.random((n_data, n_vbin))
        ref = reference(y, f)
        for chunk in (1, 7, n_orb, n_orb + 100):   # ragged, exact, oversized
            got = chunked(y, f, chunk)
            assert got.shape == ref.shape, (got.shape, ref.shape)
            assert np.allclose(got, ref), f'mismatch at {chunk=}'
    print("chunked LOSVD rebin == einsum reference, OK")


if __name__ == "__main__":
    demo()
