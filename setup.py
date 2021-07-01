#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
from __future__ import absolute_import, print_function

import os
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup

if 'TOX_ENV_NAME' not in os.environ \
        or 'TOX_ENV_NAME' in os.environ and 'TOXTESTING' in os.environ:

    from pybind11.setup_helpers import Pybind11Extension, build_ext

    ext_modules = [
        Pybind11Extension(
            "idpcpp",
            #sorted(glob(join('src', 'idpconfgen', 'cpp', '*.cpp'))),
            [join('src', 'idpconfgen', 'cpp', 'idpcpp.cpp')],
            )
        ]
    cmdclass={'build_ext': build_ext}

else:
    ext_modules = []
    cmdclass = {}


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()


long_description = '{}\n{}'.format(
    re.compile(
        '^.. start-badges.*^.. end-badges',
        re.M | re.S,
        ).sub(
            '',
            read('README.rst'),
            ),
    re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read(join('CHANGELOG.rst')))
    )

setup(
    name='idpconfgen',
    version='0.0.22',
    description='Generates IDP conformers.',
    long_description=long_description,
    author='Julie Forman-Kay Lab',
    author_email='forman@sickkids.ca',
    url='https://github.com/julie-forman-kay-lab/IDPConformerGenerator',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    project_urls={
        'webpage': 'https://github.com/julie-forman-kay-lab/IDPConformerGenerator',
        # 'Documentation': 'https://python-project-skeleton.readthedocs.io/',
        # 'Changelog': 'https://python-project-skeleton.readthedocs.io/en/latest/changelog.html',
        # 'Issue Tracker': 'https://github.com/joaomcteixeira/python-project-skeleton/issues',
        },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
        ],
    python_requires='>=3.7.*,<4',
    install_requires=[
        #'libfuncpy>=0.0.3',
        #'numpy>=1,<2',
        #'numba>=0.53.0',
        #'scipy>=1,<2',
        #'pybind11>=2,<3',
        ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
        },
    setup_requires=[
        #   'pytest-runner',
        #   'setuptools_scm>=3.3.1',
        ],
    entry_points={
        'console_scripts': [
            'idpconfgen = idpconfgen.cli:maincli',
            ]
        },
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    )
