import os, sys
from setuptools import setup, find_packages

with open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='BIO-PEPPA',
    version='1.1.6',
    #scripts=['PEPPA.py'] ,
    author="Zhemin Zhou",
    author_email="zhemin.zhou@warwick.ac.uk",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zheminzhou/PEPPA",
    packages = ['PEPPA'],
    package_dir = {'PEPPA':'.'},
    keywords=['bioinformatics', 'microbial', 'genomics', 'MLST', 'pan-genome'],
    install_requires=['ete3>=3.1.1', 'numba>=0.47.0', 'numpy>=1.18.1', 'pandas>=1.0.1', 'scipy>=1.3.2'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'PEPPA = PEPPA.PEPPA:ortho',
            'PEPPA_parser = PEPPA.PEPPA_parser:PEPPA_parser',
    ]},
    package_data={'PEPPA': ['LICENSE', 'README.*', 'dependencies/*', 'modules/*.py']},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent",
    ],
 )

