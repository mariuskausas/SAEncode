"""
SAEncode
Python implementation of Structural Alphabet encoding of molecular dynamics trajectories.
"""

from setuptools import setup
from setuptools import setup, find_packages

setup(
    # Self-descriptive entries which should always be present
    name='saencode',
    author='Marius Kausas',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    package_data={'saencode': ['alphabet/*.pdb']},
    version=0.1,
)
