from pathlib import Path

from setuptools import setup


def run_setup(version):
    cdir = Path(__file__).parent.absolute()
    with open(cdir.joinpath("README.md")) as f:
        long_description = f.read()

    with open(cdir.joinpath("requirements.txt")) as f:
        requirements = [line.strip() for line in f.readlines()]

    # We need to define all the modules explicitly for setup()
    modules = [
        "cryptography318.dlp",
        "cryptography318.factor",
        "cryptography318.linalg",
        "cryptography318.prime",
        "cryptography318.utils"
    ]

    setup(
        script_name='setup.py',
        script_args=['sdist', 'bdist_wheel'],
        name='cryptography318',
        long_description=long_description,
        long_description_content_type='text/markdown',
        version=version,
        description='A set of functions useful in mathematical cryptography',
        url='https://github.com/aarpyy/Cryptography',
        author='Andrew Carpenter',
        author_email='andrewcarp00@gmail.com',
        project_urls={
            'Source': 'https://github.com/aarpyy/Cryptography'
        },
        packages=["cryptography318"] + modules,
        install_requires=requirements,
        python_requires=">=3.10",  # Python 3.10 now required since zip(*, strict=True) is used
        classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Education',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
        ],
    )
