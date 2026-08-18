"""Deliberately break the bulk parser and confirm the test suite notices.

A parser test that cannot fail is worthless, so mutate each load-bearing line
and check the suite goes red.
"""
import subprocess
import sys

P = 'dynamite/orblib.py'
MUTS = [
    ("walker stride 6->4",
     "            p += 6\n",
     "            p += 4\n"),
    ("restore the bogus header skip",
     "        p = 0\n        start = [0] * n_pairs\n        ivmin0",
     "        p = (4 + mv[0] + 4) // 4\n        start = [0] * n_pairs\n        ivmin0"),
    ("payload decoded as C order not F",
     "            i0 = within % n0_e\n            i1 = within // n0_e",
     "            i1 = within % n0_e\n            i0 = within // n0_e"),
    ("swap row/col destination",
     "            row = np.repeat(ivmin0[k] + centre0[kin], nvs) + i0\n"
     "            col = np.repeat(ivmin1[k] + centre1[kin], nvs) + i1",
     "            row = np.repeat(ivmin0[k] + centre0[kin], nvs) + i1\n"
     "            col = np.repeat(ivmin1[k] + centre1[kin], nvs) + i0"),
    ("forget idx_ap_reset in the pm scatter",
     "            ap_e = np.repeat(ap - idx_ap_reset[kin], nvs)",
     "            ap_e = np.repeat(ap, nvs)"),
    ("wrong fast-axis length n0",
     "                n0[k] = i_max0 - i_min0 + 1",
     "                n0[k] = i_max1 - i_min1 + 1"),
    ("1d parser ignores the aperture subset",
     "            ap = ap_global[k - orb * n_ap_file]\n            kin = kin_idx_per_ap[ap]\n"
     "            # expand each pair into its nv individual values",
     "            ap = k - orb * n_ap_file\n            kin = kin_idx_per_ap[ap]\n"
     "            # expand each pair into its nv individual values"),
]

orig = open(P).read()
results = []
try:
    for name, old, new in MUTS:
        if orig.count(old) != 1:
            results.append(("SKIP", name, f"anchor appears {orig.count(old)}x"))
            continue
        open(P, 'w').write(orig.replace(old, new))
        r = subprocess.run([sys.executable, 'dev_tests/test_hist_bulk_read.py'],
                           capture_output=True, text=True)
        results.append(("CAUGHT" if r.returncode else "MISSED", name, ""))
finally:
    open(P, 'w').write(orig)

for status, name, extra in results:
    print(f"{status:7s} {name} {extra}")
print("\n(file restored)")
