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

    # Make sure we have setuptools and twine before running setup.py
    for pkg in ("setuptools", "twine", "wheel"):
        if importlib.util.find_spec(pkg) is None:
            sp.run([sys.executable, "-m", "pip", "install", pkg], cwd=cdir)

    # Run setup.py
    sp.run([sys.executable, "setup.py", "sdist", "bdist_wheel"], cwd=cdir)

    # If test, upload to test.pypi.org
    if len(sys.argv) > 1 and sys.argv[1] in ("-t", "--test"):
        print("Uploading package to test.pypi.org")
        sp.run(["twine", "upload", "--repository-url", "https://test.pypi.org/legacy/", "dist/*"], cwd=cdir)
    else:
        print("Uploading package to pypi.org")
        # Otherwise, upload to pypi.org
        sp.run(["twine", "upload", "dist/*"], cwd=cdir)


if __name__ == "__main__":
    main()
