from setuptools import setup
from os import path


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='cryptography318',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.1.5',
    description='A set of functions useful in cryptography.',
    url='https://github.com/aarpyy/Cryptography',
    author='Andrew Carpenter',
    author_email='acarpent@oberlin.edu',
    packages=['cryptography318'],
    install_requires=['numpy'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
