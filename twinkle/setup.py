from setuptools import setup, find_packages
import os

# def read(fname):
#     return open(os.path.join(os.path.dirname(__file__), fname)).read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='twinkle',
    version='1.0.0',
    description='A custom package for calculating stellar spectral energy distributions (SED) using broadband fluxes and grid models. Also has the capability to calculate excess flux.',
    author='Rahul I. Patel',
    url='https://github.com/astropatel/twinkle',
    author_email='ri.patel272@gmail.com',
    include_package_data=True,
    packages=find_packages(where='twinkle'),
    package_dir={"":"twinkle"},
    install_requires=required,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    python_requires='>3.8'
)
