"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-chem-iedb',
    version='3.0.0',
    description='Scipion framework plugin for the use of iedb for T-cell epitope prediction',
    long_description=long_description,
    url='https://github.com/scipion-chem/scipion-chem-tepitool',
    author='Daniel Del Hoyo Gomez',
    author_email='scipion@cnb.csic.es',
    keywords='scipion epitope t-cell',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    package_data={
       'iedb': ['iedb_logo.png'],
    },
    entry_points={
        'pyworkflow.plugin': 'iedb = iedb'
    }
)
