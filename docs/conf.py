# -*- coding: utf-8 -*-
"""Configuration file for Sphinx."""
from __future__ import unicode_literals

import os
import sys

import mock


mock_modules = [
    'numba',
    'scipy',
    'scipy.spatial',
    ]

for modulename in mock_modules:
    sys.modules[modulename] = mock.Mock()

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    ]

todo_include_todos = True

exclude_patterns = [
    'i_*',
    ]

if os.getenv('SPELLCHECK'):
    extensions += 'sphinxcontrib.spelling',
    spelling_show_suggestions = True
    spelling_lang = 'en_US'
    # https://sphinxcontrib-spelling.readthedocs.io/en/latest/customize.html
    spelling_word_list_filename = ['spelling_wordlist.txt']

source_suffix = '.rst'
master_doc = 'index'
project = 'IDPConformerGenerator'
year = '2022'
author = 'Julie Forman-Kay Lab'
copyright = '{0}, {1}'.format(year, author)
version = release = '0.7.31'

pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/julie-forman-kay-lab/IDPConformerGenerator/issues/%s', '#'),  # noqa: E501
    'pr': ('https://github.com/julie-forman-kay-lab/IDPConformerGenerator/pulls/%s', 'PR #'),  # noqa: E501
    }


# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
    '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
    }
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
