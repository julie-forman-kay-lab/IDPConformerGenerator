"""
IDP Conformer Generator.

Generates conformers for Intrinsically Disordered Proteins.
"""
import argparse
import logging
import os
from pathlib import Path as _Path
import sys


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter('%(message)s')

log.addHandler(_ch)


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


class Path(type(_Path())):
    """
    Define a common Path to string interface.

    Avoids using os.fspath around libs.
    """
    def str(self):
        """Return string version of Path."""
        return os.fspath(self)


__version__ = '0.0.4'
