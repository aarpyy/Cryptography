from setuptools import setup

setup(
    name='cryptography318',
    version='0.1.0',
    description='A set of functions useful in cryptography',
    url='https://github.com/aarpyy/Cryptography',
    author='Andrew Carpenter',
    author_email='acarpent@oberlin.edu',
    license='BSD 2-clause',
    packages=['cryptography318'],
    install_requires=['mpi4py>=2.0',
                      'numpy',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
