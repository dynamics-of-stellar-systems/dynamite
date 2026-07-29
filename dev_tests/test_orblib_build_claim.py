"""Standalone check for LegacyOrbitLibrary's orbit-library build lock."""

import os
import tempfile
from types import MethodType


def main():
    import dynamite.orblib

    cls = dynamite.orblib.LegacyOrbitLibrary
    logger = type("L", (), {"warning": lambda s, m: None})()

    with tempfile.TemporaryDirectory() as tmp:
        mod_dir = tmp + "/"

        # Build a fake that delegates all three lock methods to the class
        class Fake:
            pass

        Fake.claim_orblib_build = lambda s: MethodType(cls.claim_orblib_build, s)()
        Fake.release_orblib_build_claim = lambda s: MethodType(cls.release_orblib_build_claim, s)()
        Fake._orblib_lock_is_stale = lambda s, lp: MethodType(cls._orblib_lock_is_stale, s)(lp)

        obj = Fake()
        obj.mod_dir = mod_dir
        obj.logger = logger

        assert obj.claim_orblib_build() is True, "first claim should succeed"
        assert obj.claim_orblib_build() is False, "second concurrent claim should be refused"

        obj.release_orblib_build_claim()
        assert obj.claim_orblib_build() is True, "claim after release should succeed"

        # stale lock detection
        lock_path = mod_dir + "datfil/building.lock"
        os.makedirs(mod_dir + "datfil", exist_ok=True)
        dead_pid = 2**30
        with open(lock_path, "w") as f:
            f.write(str(dead_pid))
        assert obj._orblib_lock_is_stale(lock_path) is True
        assert obj.claim_orblib_build() is True, "stale lock should be retaken"
        obj.release_orblib_build_claim()

        print("test_orblib_build_claim: all checks passed")


if __name__ == "__main__":
    main()
