from setuptools import setup
from pathlib import Path


cdir = Path(__file__).parent.absolute()
with open(cdir.joinpath('README.md')) as f:
    long_description = f.read()


setup(
    name='cryptography318',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.2.1',
    description='A set of functions useful in cryptography and linear algebra',
    url='https://github.com/aarpyy/Cryptography',
    author='Andrew Carpenter',
    author_email='acarpent@oberlin.edu',
    packages=['cryptography318'],
    install_requires=['numpy',
                      'sympy'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Education',
        'License :: Public Domain',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows :: Windows 10',
        'Programming Language :: Python :: 3.10'        # Python 3.10 now required since zip(*, strict=True) is used
    ],
)
