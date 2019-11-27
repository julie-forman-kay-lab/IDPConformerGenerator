"""Library for the Conformer builder."""
import numpy as np

from idpconfgen.core import definitions as DEFS


class ConformerTemplate:
    """
    A conformer template for building.

    Used by the :class:`ConformerBuilder` as a molecule coordinate
        system.
    """
    def __init__(self, seq):
        
        self._seq = self._parse_seq(seq)
        self._coords = np.full(
            len(self.seq) * len(DEFS.backbone_atoms) + 1, # the oxygen from COO
            np.nan,
            )

    @property
    def seq(self):
        return self._seq

    @property
    def coords(self):
        return self._coords

    @staticmethod
    def _parse_seq(seq):
        """
        seq (list or str) -> tuple of 3 char strings
        """
        try:
            sequence = tuple(DEFS.aa1to3[c] for c in seq)
        except KeyError:
            # the sequece is not a sequence of 1 letter codes
            if all(c in DEFS.aa3to1.keys() for c in seq):
                sequence = tuple(seq)
            else:
                raise ValueError('seq aminoacids codes are not valid')
        return sequence
