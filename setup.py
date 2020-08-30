#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

# Fetch the version
exec(open('workflowhub/version.py').read())

setup(
    name='workflowhub',
    version=str(__version__),
    license='GPLv3',
    author='WorkflowHub team',
    author_email='support@workflowhub.org',
    description='Community Framework for Enabling Scientific Workflow Research and Education',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/workflowhub/workflowhub',
    packages=find_packages(),
    install_requires=[
        'jsonschema',
        'matplotlib',
        'networkx',
        'numpy',
        'pygraphviz',
        'python-dateutil',
        'requests',
        'scipy',
        'setuptools'
    ],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Documentation :: Sphinx',
        'Topic :: System :: Distributed Computing'
    ],
    python_requires='>=3.5',
)
