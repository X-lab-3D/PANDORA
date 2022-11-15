# -*- coding: utf-8 -*-
import os
from pathlib import Path
from os.path import exists
import json

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit PANDORA/__version__.py
version = {}
with open(os.path.join(here, 'PANDORA', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.md') as readme_file:
    readme = readme_file.read()

def post_install_jobs():
    user_folder_path = Path(__file__).parents[0]

    if exists('config.json'):
        with open('config.json') as f:
            data = json.load(f)
            data_folder = data['data_folder_name']
    else:
        data_folder = 'default'

    dirs = [
            f'{user_folder_path}/Databases', 
            f'{user_folder_path}/Databases/{data_folder}',
            f'{user_folder_path}/Databases/{data_folder}/mhcseqs', 
            f'{user_folder_path}/Databases/{data_folder}/BLAST_databases',
            f'{user_folder_path}/Databases/{data_folder}/PDBs',
            f'{user_folder_path}/Databases/{data_folder}/PDBs/pMHCI', 
            f'{user_folder_path}/Databases/{data_folder}/PDBs/pMHCII',
            f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad', 
            f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad/pMHCI',
            f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad/pMHCII', 
            f'{user_folder_path}/Databases/{data_folder}/PDBs/IMGT_retrieved',
            f'{user_folder_path}/Databases/{data_folder}/outputs',
            f'{user_folder_path}/test/',
            f'{user_folder_path}/test/test_data',
            f'{user_folder_path}/test/test_data/PDBs/Bad',
            f'{user_folder_path}/test/test_data/PDBs/Bad/pMHCI',
            f'{user_folder_path}/test/test_data/PDBs/Bad/pMHCII', 
            ]

    for D in dirs:
        try:
            os.mkdir(D)
        except OSError as e:
            print(f'Could not make directory: {D} \n Reason: {e}')

    try:
        print('Downloading pre-built database from zenodo...')
        os.popen(f'wget https://sandbox.zenodo.org/record/1129456/files/default.tar.gz?download=1 -O {user_folder_path}/Databases/default.tar.gz').read()
        print('Copying the database')
        os.popen(f'tar -xzvf {user_folder_path}/Databases/default.tar.gz -C {user_folder_path}/Databases/{data_folder}').read()
        os.popen(f'rm {user_folder_path}/Databases/default.tar.gz').read()
    except Exception as e:
        print(f'WARNING: received error while installing database: {e}')
        print('To be able to use PANDORA you will have to generate a new database. Please follow the instructions in the README.')

class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)
        post_install_jobs()

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        post_install_jobs()


setup(
    name='PANDORA',
    version=version['__version__'],
    description='Model peptide-MHC complexes using anchor distance restrains in MODELLER',
    long_description=readme + '\n\n',
    long_description_content_type='text/markdown',
    author='Dario Marzella, Farzaneh Parizi, Li Xue',
    url='https://github.com/X-lab-3D/PANDORA/tree/master',
    project_urls={
        'Source Code': 'https://github.com/X-lab-3D/PANDORA/tree/master',
        'Issue tracker': 'https://github.com/X-lab-3D/PANDORA/issues'
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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],

    install_requires=[
       'Biopython',
       'pdb2sql',
       'joblib'
       ],

    extras_require={
        'doc': ['recommonmark', 'sphinx', 'sphinx_rtd_theme'],
        'test': ['pytest', 'pytest-runner', 'pytest-cov',
                 'coverage', 'coveralls', 'pycodestyle']
    },
    
    cmdclass={
        'develop': PostInstallCommand,
        'install': PostInstallCommand,
    },
)