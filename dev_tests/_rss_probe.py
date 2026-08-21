"""Heap-reclaim probe + RSS sampler for the weight-solve memory work.

Two independent tools:

- ``reclaim_probe``: allocate arrays of given sizes, free them, report how
  much VmRSS the OS took back after del + gc.collect() (+ malloc_trim(0)).
  Answers whether Linux glibc reclaims huge numpy allocations here; the
  earlier macOS investigation measured ~0% reclaim at every size.
- ``RssSampler``: daemon thread sampling VmRSS/VmHWM (/proc/self/status) and
  Anonymous/PSS (/proc/self/smaps_rollup) into a CSV tagged with the current
  phase name.

Run standalone:
    python dev_tests/_rss_probe.py --sizes 1 8 64   # GiB
    python dev_tests/_rss_probe.py --selfcheck
"""

import argparse
import ctypes
import gc
import os
import threading
import time

import numpy as np


def _read_status():
    out = {}
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith(("VmRSS:", "VmHWM:")):
                key, val = line.split(":", 1)
                out[key] = int(val.strip().split()[0])  # kB
    return out


def _read_smaps():
    out = {}
    try:
        with open("/proc/self/smaps_rollup") as fh:
            for line in fh:
                if line.startswith(("Anonymous:", "Pss:")):
                    key, val = line.split(":", 1)
                    out[key] = int(val.strip().split()[0])  # kB
    except OSError:
        pass
    return out


class RssSampler(threading.Thread):
    """Sample this process's RSS into a CSV until the context exits.

    The ``phase`` attribute is mutable from the sampling loop's own process;
    each row is tagged with whatever it holds at sample time."""

    def __init__(self, path, interval=1.0):
        super().__init__(daemon=True)
        self.path = path
        self.interval = interval
        self.phase = "init"
        self._stop_evt = threading.Event()

    def __enter__(self):
        self._fh = open(self.path, "w", buffering=1)
        self._fh.write("epoch_s,phase,vmrss_kb,vmhwm_kb,anon_kb,pss_kb\n")
        self._row()  # guarantee one row at the entering phase
        self.start()
        return self

    def __exit__(self, *exc):
        self._stop_evt.set()
        self.join(timeout=2 * self.interval + 1)
        self._fh.close()

    def _row(self):
        st = _read_status()
        sm = _read_smaps()
        self._fh.write(
            f"{time.monotonic():.1f},{self.phase},"
            f"{st.get('VmRSS', -1)},{st.get('VmHWM', -1)},"
            f"{sm.get('Anonymous', -1)},{sm.get('Pss', -1)}\n"
        )

    def run(self):
        while not self._stop_evt.is_set():
            self._row()
            self._stop_evt.wait(self.interval)


def _malloc_trim():
    try:
        return bool(ctypes.CDLL("libc.so.6").malloc_trim(0))
    except OSError:
        return False


def reclaim_probe(sizes_gib):
    """Allocate/free arrays of the given sizes and report RSS at each stage."""
    gib = 2**30
    print(
        f"{'GiB':>6} {'rss@alloc':>11} {'after_del':>11} "
        f"{'after_gc':>11} {'after_trim':>11} reclaimed%"
    )
    for s in sizes_gib:
        base = _read_status()["VmRSS"]
        arr = np.empty(int(s * gib // 8), dtype=np.float64)
        arr[:] = 1.0  # fill: force every page resident (bandwidth-bound)
        rss_alloc = _read_status()["VmRSS"]
        del arr
        rss_del = _read_status()["VmRSS"]
        gc.collect()
        rss_gc = _read_status()["VmRSS"]
        _malloc_trim()
        rss_trim = _read_status()["VmRSS"]
        gained = rss_alloc - rss_trim  # kB handed back
        pct = 100.0 * gained / max(rss_alloc - base, 1)
        print(
            f"{s:>6} {rss_alloc / 2**20:>10.2f}G {rss_del / 2**20:>10.2f}G "
            f"{rss_gc / 2**20:>10.2f}G {rss_trim / 2**20:>10.2f}G {pct:>8.1f}%"
        )


def _selfcheck():
    with RssSampler("/tmp/_rss_selfcheck.csv", interval=0.2) as samp:
        samp.phase = "alloc"
        junk = np.empty(2**28, dtype=np.uint8)
        junk[:] = 1
        time.sleep(0.6)
        samp.phase = "freed"
        del junk
        gc.collect()
        time.sleep(0.6)
    lines = open("/tmp/_rss_selfcheck.csv").read().strip().splitlines()
    assert len(lines) > 4, lines
    phases = {ln.split(",")[1] for ln in lines[1:]}
    assert {"init", "alloc", "freed"} <= phases, phases
    print(f"selfcheck OK: {len(lines) - 1} rows, phases={sorted(phases)}")
    os.remove("/tmp/_rss_selfcheck.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=float, nargs="*", default=[1, 8, 64])
    ap.add_argument("--selfcheck", action="store_true")
    args = ap.parse_args()
    if args.selfcheck:
        _selfcheck()
    else:
        reclaim_probe(args.sizes)
