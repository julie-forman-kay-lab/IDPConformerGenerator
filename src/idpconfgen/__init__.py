"""
IDP Conformer Generator.

Generates conformers for Intrinsically Disordered Proteins.
"""

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter('%(message)s')

log.addHandler(_ch)


__version__ = '0.0.4'
