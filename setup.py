# -*- coding: utf-8 -*-
import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit PANDORA/__version__.py
version = {}
with open(os.path.join(here, 'PANDORA', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='PANDORA',
    version=version['__version__'],
    description='Model peptide-MHC I complexes using anchor distance restrains in MODELLER',
    long_description=readme + '\n\n',
    long_description_content_type='text/markdown',
    author='Farzaneh Meimandi Parizi, Dario Marzella, Li Xue',
    url='https://github.com/DarioMarzella/pMHC_Modelling/tree/master',
    project_urls={
        'Source Code': 'https://github.com/DarioMarzella/pMHC_Modelling/tree/master',
        'Issue tracker': 'https://github.com/DarioMarzella/pMHC_Modelling/issues'
    },
    packages=find_packages(),
    include_package_data=True,
    license="Apache Software License 2.0",
    keywords='PANDORA',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    install_requires=[
       'Biopython',
       'pdb_tools',
       'pdb2sql',
       'muscle',
       'matplotlib',
       'dill',
       'pathos'
       ],

    extras_require={
        'doc': ['recommonmark', 'sphinx', 'sphinx_rtd_theme'],
        'test': ['pytest', 'pytest-runner',
                 'coverage', 'coveralls', 'pycodestyle']
    }
)
