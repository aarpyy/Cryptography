from pathlib import Path
import subprocess as sp
import importlib.util
import sys


def main():
    cdir = Path(__file__).parent

    # Remove all previous package build stuff
    for f in cdir.iterdir():
        if f.name.endswith(".egg-info") or (f.name in ("dist", "build") and f.is_dir()):
            sp.run(["rm", "-rf", f.absolute().name])

    # Make sure we have setup tools before running setup.py
    if importlib.util.find_spec("setuptools") is None:
        sp.run([sys.executable, "-m", "pip", "install", "setuptools"], cwd=cdir)
    if importlib.util.find_spec("twine") is None:
        sp.run([sys.executable, "-m", "pip", "install", "twine"], cwd=cdir)

    # Run setup.py
    sp.run(["python3", "setup.py", "sdist", "bdist_wheel"], cwd=cdir)

    if len(sys.argv) > 1 and sys.argv[1] in ("-t", "--test"):
        # And upload to testpypi
        sp.run(["twine", "upload", "--repository-url", "https://test.pypi.org/legacy/", "dist/*"], cwd=cdir)
    else:
        # And upload to pypi
        sp.run(["twine", "upload", "dist/*"], cwd=cdir)


if __name__ == "__main__":
    main()
