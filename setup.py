#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', ]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Adam James Finley",
    author_email='af472@exeter.ac.uk',
    python_requires='!=2.7, >=3.0., !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="PLUTO stellar wind simulation visualisation tools and analysis scripts.",
    entry_points={
        'console_scripts': [
            'pluto_stellarwindpy=pluto_stellarwindpy.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pluto_stellarwindpy',
    name='pluto_stellarwindpy',
    packages=find_packages(include=['pluto_stellarwindpy', 'pluto_stellarwindpy.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/af472/pluto_stellarwindpy',
    version='0.1.0',
    zip_safe=False,
)
