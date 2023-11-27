import os
from setuptools import setup


def get_version():
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sphae', 'sphae.VERSION')) as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name='sphae',
    description="Assembling pure culture phages from both Illumina and Nanopore sequencing technology",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Bhavya Papudeshi",
    author_email="npbhavya13@gmail.com",
    data_files=get_data_files(),
    py_modules=['sphae'],
    install_requires=[
        "snaketool-utils>=0.0.4",
        "snakemake>=7.14.0",
        "pyyaml>=6.0",
        "Click==8.1.3",
        "metasnek>=0.0.4",
        "attrmap>=0.0.7",
        "biopython>=1.8.1",
        "pandas"
    ],
    entry_points={
        'console_scripts': [
            'sphae=sphae.__main__:main'
        ]
    },
    include_package_data=True,
)
