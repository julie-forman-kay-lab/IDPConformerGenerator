"""Test Structure class."""
from pathlib import Path

import numpy as np
import pytest
import hypothesis.strategies as st
from hypothesis import given

from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libstructure import (
    gen_empty_structure_data_array,
    parse_pdb_to_array,
    parse_cif_to_array,
    Structure,
    )

from . tcommons import (
    pdb_example,
    cif_example,
    cif_noatomsite,
    cif_nohash,
    cif_EOF,
    pdb_saved,
    )


@pytest.mark.parametrize(
    'is_func,structure,expected',
    [
        (is_cif, cif_example, True),
        (is_cif, pdb_example, False),
        (is_pdb, pdb_example, True),
        ]
    )
def test_is_function(is_func, structure, expected):
    """
    """
    assert is_func(structure.read_text()) == expected


@given(st.integers(min_value=1, max_value=100_000))
def test_gen_array_data(number):
    """
    """
    result = gen_empty_structure_data_array(number)
    assert result.shape == (number, len(PDBParams.atom_slicers))
    assert result.dtype == '<U8'


def test_parse_pdb_to_array():
    """
    """
    data_array = parse_pdb_to_array(pdb_example.read_text())
    assert data_array.shape == (5385, len(PDBParams.atom_slicers))
    assert data_array.dtype == np.dtype('<U8')


def test_parse_pdb_to_array_ValueError():
    """."""
    with pytest.raises(ValueError):
        parse_pdb_to_array(pdb_example.read_text(), which='none')


# TODO
def test_parse_cif_to_array():
    """
    """
    data_array = parse_cif_to_array(cif_example.read_text())
    assert data_array.shape == (5385, len(PDBParams.atom_slicers))
    assert data_array.dtype == np.dtype('<U8')


@pytest.mark.parametrize(
    'ciffile',
    [
        cif_noatomsite,
        cif_nohash,
        cif_EOF,
        ]
    )
def test_parse_cif_to_array_noatoms_error(ciffile):
    """/"""
    with pytest.raises(EXCPTS.CIFFileInvalidError):
        parse_cif_to_array(ciffile.read_text())


@pytest.fixture()
def CIFParser_fixture(request):
    return CIFParser(cif_example.read_text())


@pytest.mark.parametrize(
    'example,parser',
    [
        (pdb_example, parse_pdb_to_array),
        (cif_example, parse_cif_to_array),
        ]
    )
def test_Structure_dispatch(example, parser):
    """
    Test Structure dispatch.

    Structure should return the correct subclass.
    """
    result = Structure._detect_structure_type(example.read_text())
    assert result == parser


@pytest.fixture(params=[
    pdb_example,
    cif_example,
    ])
def Structure_built(request):
    s = Structure(request.param.read_text())
    s.build()
    return s


@pytest.mark.parametrize(
    'func,number_of_lines',
    [
        (lambda x: int(x[PDBParams.acol.resseq]) < 10, 96),
        ]
    )
def test_Structure_1(Structure_built, func, number_of_lines):
    """
    """
    s = Structure_built
    s.add_filter(func)
    assert len(list(s.filtered_atoms)) == number_of_lines


@pytest.mark.parametrize(
    'record',
    [
        ('ATOM', 'HETATM'),
        'ATOM',
        'HETATM',
        ]
    )
def test_Structure_filter_record_name(Structure_built, record):
    """
    """
    Structure_built.add_filter_record_name(record)
    lines = list(Structure_built.filtered_atoms)
    assert all(l[0].strip() in record for l in lines)


@pytest.mark.parametrize(
    'chain, expected',
    [
        ('A',2693),
        ]
    )
def test_Structure_filter_chain(Structure_built, chain, expected):
    """
    """
    Structure_built.add_filter_chain(chain)
    assert len(list(Structure_built.filtered_atoms)) == expected


def test_Structure_chain_set(Structure_built):
    """
    """
    assert Structure_built.chain_set == set(['A', 'B'])


def test_Structure_save(Structure_built):
    """Test save functionality."""
    Structure_built.add_filter(lambda x: int(x[PDBParams.acol.resseq]) < 11)
    Structure_built.add_filter_chain('A')

    fout = Path('pdb_testing.pdb')
    Structure_built.write_PDB(fout)
    result = fout.read_text()
    expected = pdb_saved.read_text()
    fout.unlink()
    assert result == expected


def test_Structure_save_empty_filter_error(Structure_built):
    """
    """
    Structure_built.add_filter_chain('Z')
    with pytest.raises(EXCPTS.EmptyFilterError):
        Structure_built._make_pdb()
