"""Test Structure class."""
import inspect
from pathlib import Path

import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import given

from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libpdb
from idpconfgen.libs import libstructure

from . tcommons import (
    cif_example,
    pdb_example,
    pdb_models,
    pdb_models_output,
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


@pytest.fixture(
    params=[
        (pdb_example, libstructure.parse_pdb_to_array),
        (cif_example, libstructure.parse_cif_to_array),
        ],
    )
def parsed_array_examples(request):
    """Parse functions."""
    return request.param[1](request.param[0].read_text())


def test_parsed_array_examples_shape(parsed_array_examples):
    """Test shape of the parsed array."""
    assert parsed_array_examples.shape == (545, len(libpdb.atom_slicers))


def test_parsed_array_examples_dtype(parsed_array_examples):
    """Test array data is U8."""
    assert parsed_array_examples.dtype == np.dtype('<U8')


def test_parse_pdb_to_array_ValueError():
    """Test parse_pdb_to_array raises ValueError."""
    with pytest.raises(ValueError):
        libstructure.parse_pdb_to_array(pdb_example.read_text(), which='none')


@pytest.fixture(params=[
    cif_example,
    pdb_example,
    ])
def fix_Structure(request):
    """Instantiate Structure object from file."""
    s = libstructure.Structure(request.param.read_text())
    return s


@pytest.fixture(params=[
    cif_example,
    pdb_example,
    ])
def fix_Structure_build(request):
    """Instantiate Structure object from file builded."""
    s = libstructure.Structure(request.param.read_text())
    s.build()
    return s


def test_structure_parsers(fix_Structure):
    """Test has structure_parser function."""
    assert inspect.isfunction(fix_Structure._structure_parser)


def test_has_datastr(fix_Structure):
    """Test has attribute."""
    assert fix_Structure._datastr


def test_not_data_array(fix_Structure):
    """Test AttributeError."""
    with pytest.raises(EXCPTS.NotBuiltError):
        fix_Structure.data_array


def test_Structure_build_1(fix_Structure_build):
    """Build method."""
    with pytest.raises(AttributeError):
        fix_Structure_build._datastr


def test_Structure_build_2(fix_Structure_build):
    """Test data array is numpy array."""
    assert isinstance(fix_Structure_build.data_array, np.ndarray)


def test_Structure_fasta_A(fix_Structure_build):
    """Test FASTA property."""
    fix_Structure_build.add_filter_chain('A')
    fix_Structure_build.add_filter_record_name('ATOM')
    fasta = fix_Structure_build.fasta
    assert fasta == {'A': 'AYIAKQRQISFVKSHF'}
    # this functionality is used in cli_fastaext
    assert next(iter(fasta.values())) == 'AYIAKQRQISFVKSHF'


def test_Structure_fasta_CIF_all():
    """Test FASTA property."""
    s = libstructure.Structure(cif_example)
    s.build()
    expected = {
        'A': 'AYIAKQRQISFVKSHF',
        'B': 'AYIAKQRQISFVKSHFS',
        'C': 'N',
        'D': 'X',
        'E': 'N',
        'F': 'X',
        'G': 'X',
        'H': 'X',
        }
    fasta = s.fasta
    assert fasta == expected
    assert next(iter(fasta.values())) == 'AYIAKQRQISFVKSHF'


def test_Structure_build_3(fix_Structure_build):
    """Test shape of data_array."""
    expected = (545, len(libpdb.atom_slicers))
    assert fix_Structure_build.data_array.shape == expected


@pytest.fixture(
    params=[
        (True, 545),
        (False, 0),
        ])
def bool_(request):
    """Boolean filters and corresponding structure number of atoms."""
    return request.param


@pytest.fixture
def filtered(request, fix_Structure_build, bool_):
    """Test Structure build filter false."""
    fix_Structure_build.add_filter(lambda x: bool_[0])
    return fix_Structure_build, bool_[1]


def test_false_num_of_filters(filtered):
    """Test number of filters."""
    assert len(filtered[0].filters) == 1


def test_false_len_filtered_atoms(filtered):
    """Test len of filtered atoms."""
    assert len(list(filtered[0].filtered_atoms)) == filtered[1]


def test_false_pop_last_filter(filtered):
    """Test pop_last_filter."""
    filtered[0].pop_last_filter()
    assert len(filtered[0].filters) == 0


@pytest.mark.parametrize(
    'in1,expected',
    [
        (cif_example, 142),
        (pdb_example, 275),  # :%s/\(^ATOM\|^HETATM\).\{} B \s.//n
        ]
    )
def test_Structure_filter_B(in1, expected):
    """Test Structure filter chain B."""
    s = libstructure.Structure(in1.read_text())
    s.build()
    s.add_filter_chain('B')
    result = list(s.filtered_atoms)
    assert len(result) == expected


def test_Structure_filter_ATOM(fix_Structure_build):
    """Test Structure filter ATOM."""
    fix_Structure_build.add_filter_record_name('ATOM')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 278


def test_Structure_filter_backbone(fix_Structure_build):
    """Test filter backbone atoms."""
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
    fix_Structure_build.add_filter(
        lambda x: int(x[libstructure.col_resSeq]) < 11
        )
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
