"""Setup script for the aoutools package."""
from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='aoutools',
    version='0.1.0',
    author='Jaehyun Joo',
    author_email='jaehyunjoo@outlook.com',
    description='A library of tools for analyzing All of Us data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/dokyoonkimlab/aoutools',
    packages=find_packages(),
    install_requires=[
        'hail',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
)
