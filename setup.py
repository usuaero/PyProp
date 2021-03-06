"""PyProp"""

from setuptools import setup
import os
import sys

setup(name = 'PyProp',
    version = '1.1.0',
    description = "PyProp: The ultimate electric propulsion modeling and optimization module.",
    url = 'https://github.com/usuaero/PyProp',
    author = 'usuaero',
    author_email = 'doug.hunsaker@usu.edu',
    install_requires = ['numpy>=1.18', 'scipy>=1.4', 'pytest', 'matplotlib', 'numpy-stl', 'airfoil_db @ git+https://github.com/usuaero/AirfoilDatabase.git#egg=airfoil_db-1.4.2'],
    python_requires ='>=3.6.0',
    include_package_data = True,
    license = 'MIT',
    packages = ['pyprop'],
    zip_safe = False)
