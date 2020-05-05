"""Test Structure class."""
import inspect
from pathlib import Path

import numpy as np
import pytest
import hypothesis.strategies as st
from hypothesis import given

from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libpdb
from idpconfgen.libs import libstructure

from . tcommons import (
    pdb_example,
    pdb_models,
    pdb_models_output,
    cif_example,
    cif_noatomsite,
    cif_nohash,
    cif_EOF,
    pdb_saved,
    )


@pytest.mark.parametrize(
    'data',
    [
        cif_example,
        'some string text',
        b'some bytes text',
        ]
    )
def test_get_datastr(data):
    """Test dict to get datastr."""
    assert libstructure.get_datastr(data)


def test_get_datastr_error():
    """Test dict to get datastr."""
    with pytest.raises(NotImplementedError):
        libstructure.get_datastr(int)


@given(st.integers(min_value=1, max_value=100_000))
def test_gen_empty_structure_array(number):
    """Test generation of Structure data structure."""
    result = libstructure.gen_empty_structure_data_array(number)
    assert result.shape == (number, len(libpdb.atom_slicers))
    assert result.dtype == '<U8'


@pytest.fixture(params=[
        (pdb_example, libstructure.parse_pdb_to_array),
        (cif_example, libstructure.parse_cif_to_array),
        ],
    )
def parsed_array_examples(request):
    return request.param[1](request.param[0].read_text())


def test_parsed_array_examples_shape(parsed_array_examples):
    assert parsed_array_examples.shape == (545, len(libpdb.atom_slicers))


def test_parsed_array_examples_dtype(parsed_array_examples):
    assert parsed_array_examples.dtype == np.dtype('<U8')


def test_parse_pdb_to_array_ValueError():
    """."""
    with pytest.raises(ValueError):
        libstructure.parse_pdb_to_array(pdb_example.read_text(), which='none')


@pytest.fixture(params=[
    cif_example,
    pdb_example,
    ])
def fix_Structure(request):
    s = libstructure.Structure(request.param.read_text())
    return s


@pytest.fixture(params=[
    cif_example,
    pdb_example,
    ])
def fix_Structure_build(request):
    s = libstructure.Structure(request.param.read_text())
    s.build()
    return s


def test_Structure(fix_Structure):
    """ """
    assert fix_Structure.data_array is None
    assert inspect.isfunction(fix_Structure._structure_parser)
    assert fix_Structure._datastr


def test_Structure_build_1(fix_Structure_build):
    """Build method."""
    with pytest.raises(AttributeError):
        fix_Structure._datastr


def test_Structure_build_2(fix_Structure_build):
    assert isinstance(fix_Structure_build.data_array, np.ndarray)


def test_Structure_build_3(fix_Structure_build):
    assert fix_Structure_build.data_array.shape == (545, len(libpdb.atom_slicers))


def test_Structure_build_filters_true(fix_Structure_build):
    fix_Structure_build.add_filter(lambda x: True)
    assert len(fix_Structure_build.filters) == 1
    assert len(list(fix_Structure_build.filtered_atoms)) == 545


def test_Structure_build_filters_false(fix_Structure_build):
    fix_Structure_build.add_filter(lambda x: False)
    assert len(fix_Structure_build.filters) == 1
    assert len(list(fix_Structure_build.filtered_atoms)) == 0

    fix_Structure_build.pop_last_filter()
    assert len(fix_Structure_build.filters) == 0



@pytest.mark.parametrize(
    'in1,expected',
    [
        (cif_example, 142),
        (pdb_example, 275),  # :%s/\(^ATOM\|^HETATM\).\{} B \s.//n
        ]
    )
def test_Structure_filter_B(in1, expected):
    s = libstructure.Structure(in1.read_text())
    s.build()
    s.add_filter_chain('B')
    result = list(s.filtered_atoms)
    assert len(result) == expected


def test_Structure_filter_ATOM(fix_Structure_build):
    fix_Structure_build.add_filter_record_name('ATOM')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 278


def test_Structure_filter_backbone(fix_Structure_build):
    # CIF VIM REGEX: :%s/^ATOM.\{10}\(N\|CA\|C\|O\)\s.//n
    # PDB VIM REGEX: :%s/^ATOM.\{9}\(N\|CA\|C\|O\)\s.//n
    assert libstructure.col_name == 2
    fix_Structure_build.add_filter(
        lambda x: x[libstructure.col_name] in ('N', 'CA', 'O', 'C')
        )
    fix_Structure_build.add_filter_record_name('ATOM')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 132




def test_Structure_chainset():
    """Test chain set."""
    s = libstructure.Structure(cif_example.read_text())
    s.build()
    assert s.chain_set == {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}


def test_Structure_write(fix_Structure_build):
    """Test save functionality."""
    fix_Structure_build.add_filter_record_name('ATOM')
    fix_Structure_build.add_filter(lambda x: int(x[libstructure.col_resSeq]) < 11)
    fix_Structure_build.add_filter_chain('A')

    fout = Path('pdb_testing.pdb')
    fix_Structure_build.write_PDB(fout)
    result = fout.read_text()
    expected = pdb_saved.read_text()
    fout.unlink()
    assert result == expected


def test_Structure_write_pdb_models():
    """Test Structure parsig MODELS."""
    s = libstructure.Structure(pdb_models)
    s.build()
    fout = Path('pdb_testing.pdb')
    s.write_PDB(fout)
    result = fout.read_text()
    expected = pdb_models_output.read_text()
    fout.unlink()
    assert result == expected


def test_Structure_write_empty_filter(fix_Structure_build):
    """Test save functionality."""
    fix_Structure_build.add_filter(lambda x: False)

    fout = Path('pdb_testing.pdb')
    with pytest.raises(EXCPTS.EmptyFilterError):
        fix_Structure_build.write_PDB(fout)


@pytest.mark.parametrize(
    'example,parser',
    [
        (pdb_example, libstructure.parse_pdb_to_array),
        (cif_example, libstructure.parse_cif_to_array),
        ]
    )
def test_Structure_dispatch(example, parser):
    """
    Test Structure dispatch.

    Structure should return the correct subclass.
    """
    result = libstructure.detect_structure_type(example.read_text())
    assert result == parser


def test_structure_type_error():
    """Raise error when parser not found."""
    with pytest.raises(EXCPTS.ParserNotFoundError):
        libstructure.detect_structure_type(' ')

