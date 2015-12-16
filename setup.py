#!/usr/bin/env python

from setuptools import setup

# Explicitly list bin scripts to be installed, seeing as I have a few local
# bin files that are not (yet) part of the distribution.
scripts = [
    'bin/check-pyflakes-output.sh',
    'bin/cluster-sequences-by-feature-deltas.py',
    'bin/create-database.py',
    'bin/database-details.py',
    'bin/find-matches.py',
    'bin/list-finders.py',
    'bin/list-symbols.py',
    'bin/perf.py',
    'bin/performance-test.py',
    'bin/print-distributed-database.py',
    'bin/print-sequence-features.py',
    'bin/prosite-to-json.py',
    'bin/scatter.py',
    'bin/start-distributed-database-backend.py',
    'bin/start-distributed-database.py',
    'bin/stop-distributed-database.py',
    'bin/write-hashes-string.py',
]

setup(name='light-matter',
      version='1.0.0',
      packages=['light'],
      include_package_data=True,
      url='https://github.com/acorg/light-matter',
      download_url='https://github.com/acorg/light-matter',
      author='Terry Jones, Barbara Muehlemann',
      author_email='tcj25@cam.ac.uk',
      keywords=['virus discovery'],
      classifiers=[
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      license='MIT',
      description='Sequence alignment based on predicted secondary structure',
      scripts=scripts,
      install_requires=['cffi>=1.0.0'],
      setup_requires=['cffi>=1.0.0'],
      cffi_modules=[
          './src/distance.py:ffi',
      ])
