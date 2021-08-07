from setuptools import setup
from os import path


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()


setup(
    name='cryptography318',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.1.7',
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
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
