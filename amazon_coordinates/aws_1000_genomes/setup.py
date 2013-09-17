#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="aws_1000_genomes",
    version="0.1",
    author = "James Boocock",
    author_email = "smilefreak@gmx.com",
    packages =find_packages(),
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'bam_server = aws1kg.http_bam_server:main',
            'process_bam = aws1kg.process_sample:main'
        ]
        }
    )
