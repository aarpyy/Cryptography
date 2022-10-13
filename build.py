from pathlib import Path
import setup
import argparse
import subprocess as sp
import importlib.util
import sys


def is_version_number(x: str):
    vinfo = x.split(".")
    if len(vinfo) != 3:
        return False
    return all(a.isdigit() for a in vinfo)


def main():
    parser = argparse.ArgumentParser(description="Build and upload package to PyPI")
    parser.add_argument("version", help="package version")
    parser.add_argument("-t", "--test", help="upload to test.pypi.org", action="store_true")
    args = parser.parse_args()

    if not is_version_number(args.version):
        print(f"positional version argument must be valid version number "
              f"like 'X.Y.Z', not '{args.version}'", file=sys.stderr)
        exit(2)

    cdir = Path(__file__).parent

    # Remove all previous package build stuff
    for f in cdir.iterdir():
        if f.name.endswith(".egg-info") or (f.name in ("dist", "build") and f.is_dir()):
            sp.run(["rm", "-rf", f.absolute().name])

    # Make sure we have setuptools and twine before running setup.py
    for pkg in ("setuptools", "twine", "wheel"):
        if importlib.util.find_spec(pkg) is None:
            sp.run([sys.executable, "-m", "pip", "install", pkg], cwd=cdir)

    # Run setup
    setup.run_setup(args.version)
    # sp.run([sys.executable, "setup.py", "sdist", "bdist_wheel"], cwd=cdir)

    # If test, upload to test.pypi.org
    if args.test:
        print("Uploading package to test.pypi.org")
        sp.run(["twine", "upload", "--repository-url", "https://test.pypi.org/legacy/", "dist/*"], cwd=cdir)
    else:
        print("Uploading package to pypi.org")
        # Otherwise, upload to pypi.org
        sp.run(["twine", "upload", "dist/*"], cwd=cdir)


if __name__ == "__main__":
    main()
