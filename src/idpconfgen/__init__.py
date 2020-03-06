"""
IDP Conformer Generator.

Generates conformers for Intrinsically Disordered Proteins.
"""
import logging
import os
from pathlib import Path as _Path


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter(logging.Formatter('%(message)s'))

log.addHandler(_ch)


class Path(type(_Path())):
    """
    A Path object dedicated to this software.

    Inherits from pathlib.Path.

    This creates an interface so that if new methods are required
    the Path interface does not need to be refactored across.
    """

    def str(self):
        """
        Return string version of Path.

        Avoids using os.fspath around libs.
        """
        return os.fspath(self)

    def myparents(self):
        """Return the Path to the parent folder resolved to absolute."""
        return self.resolve().parents[0]


__version__ = '0.1.0'
