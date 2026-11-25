"""Replay the 150-bound certificate against its pinned copy of upstream."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


ARCHIVE = Path(__file__).resolve().with_name("certificate.zip")
EXPECTED_SHA256 = "8387a28eb5e1e33348b8a5ebb1424d1dae8b846fa04d96b92ad7da04cfeca419"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temp-parent", type=Path, help="existing directory for the temporary extraction")
    args = parser.parse_args()

    if sys.version_info < (3, 10):
        parser.exit(1, "The certificate requires Python 3.10 or later.\n")
    if hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() != EXPECTED_SHA256:
        parser.exit(1, "certificate.zip does not match the submitted SHA-256 digest.\n")

    with tempfile.TemporaryDirectory(prefix="ps150-", dir=args.temp_parent) as folder:
        with zipfile.ZipFile(ARCHIVE) as archive:
            archive.extractall(folder)
        result = subprocess.run(
            [sys.executable, "-B", "verify_all.py", "--require-bounty"],
            cwd=folder,
        )
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
