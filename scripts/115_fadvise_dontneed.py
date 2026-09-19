#!/usr/bin/env python3
"""Best-effort, non-destructive eviction of named input files from the page cache.

This never drops host-wide caches. It applies POSIX_FADV_DONTNEED only to the regular files named
on the command line (recursing into named directories) and records every attempted file. A zero
exit status means the kernel accepted every advisory call, not that another process could not
recache a shared file before the benchmark starts.
"""

import argparse
import os
import stat
import sys
from pathlib import Path


def regular_files(paths: list[Path]):
    for path in paths:
        if path.is_dir():
            for root, dirs, files in os.walk(path):
                dirs.sort()
                for name in sorted(files):
                    candidate = Path(root, name)
                    try:
                        mode = candidate.stat().st_mode
                    except OSError as error:
                        yield candidate, error
                        continue
                    if stat.S_ISREG(mode):
                        yield candidate, None
        else:
            yield path, None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()

    advice = getattr(os, "POSIX_FADV_DONTNEED", None)
    if not hasattr(os, "posix_fadvise") or advice is None:
        print("os.posix_fadvise/POSIX_FADV_DONTNEED is unavailable", file=sys.stderr)
        return 2

    attempted = 0
    total_bytes = 0
    failures = 0
    rows: list[tuple[str, int, str]] = []
    for path, discovery_error in regular_files(args.paths):
        if discovery_error is not None:
            failures += 1
            rows.append((str(path), 0, f"stat_error:{discovery_error}"))
            continue
        try:
            size = path.stat().st_size
            fd = os.open(path, os.O_RDONLY)
            try:
                os.posix_fadvise(fd, 0, 0, advice)
            finally:
                os.close(fd)
            attempted += 1
            total_bytes += size
            rows.append((str(path), size, "accepted"))
        except OSError as error:
            failures += 1
            rows.append((str(path), 0, f"error:{error}"))

    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as handle:
        handle.write("path\tbytes\tstatus\n")
        for path, size, status in rows:
            handle.write(f"{path}\t{size}\t{status}\n")
        handle.write(f"# attempted_files\t{attempted}\n")
        handle.write(f"# advised_bytes\t{total_bytes}\n")
        handle.write(f"# failures\t{failures}\n")

    print(
        f"POSIX_FADV_DONTNEED accepted for {attempted} files / {total_bytes} bytes; "
        f"failures={failures}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
