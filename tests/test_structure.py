"""Test Structure class."""
import pytest

from idpconfgen.libs.libpdb import (
    Structure
    PDBStructure,
    CIFStructure,
    )

from . tcommons import (
    )


@pytest.mark.parametrize(
    'StructureClass,data_file',
    [
        (PDBStructure, ),
        (CIFStructure, ),
        ]
    )
def test_Structure_dispatch():
    """
    Test Structure dispatch.

    Structure should return the correct subclass.
    """


