"""Test Structure class."""
import inspect
from pathlib import Path as Path_
from types import GeneratorType

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from idpconfgen import Path
from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs import libpdb, libstructure

from .tcommons import (
    EXPL_A,
    EXPL_B,
    TmpFile,
    cif_example,
    pdb_example,
    pdb_models,
    pdb_models2,
    pdb_models_output,
    pdb_res_gap,
    pdb_saved,
    )


# # Pytest Fixtures used in this module
@pytest.fixture(
    params=[
        (pdb_example, libstructure.parse_pdb_to_array),
        (cif_example, libstructure.parse_cif_to_array),
        ],
    )
def parsed_array_examples(request):
    """Parse functions."""
    return request.param[1](request.param[0].read_text())


@pytest.fixture(
    params=[
        cif_example,
        pdb_example,
        ]
    )
def fix_Structure(request):
    """Instantiate Structure object from file."""
    s = libstructure.Structure(request.param.read_text())
    return s


@pytest.fixture
def fix_Structure_build(fix_Structure):
    """Instantiate Structure object from file builded."""
    fix_Structure.build()
    return fix_Structure


@pytest.fixture(
    params=[
        (True, 545),
        (False, 0),
        ]
    )
def bool_(request):
    """Boolean filters and corresponding expected structure number of atoms."""
    return request.param


@pytest.fixture
def filtered(request, fix_Structure_build, bool_):
    """
    Test Structure with `all` and `no` lines.

    `bool_` fixture filter selects all lines or no lines.
    """
    fix_Structure_build.add_filter(lambda x: bool_[0])
    return fix_Structure_build, bool_[1]
# end Fixture definition


# ## Start testing
# # Test module's functions
@pytest.mark.parametrize(
    'data',
    [
        cif_example,
        Path_(__file__),
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


@given(st.integers(min_value=1, max_value=100_000))
def test_gen_empty_structure_array(number):
    """Test generation of Structure data structure."""
    result = libstructure.gen_empty_structure_data_array(number)
    assert result.shape == (number, len(libpdb.atom_slicers))
    assert result.dtype == '<U8'


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
# end test module's functions


# # Test Structure __init__ definition
def test_structure_parsers(fix_Structure):
    """Test has structure_parser function."""
    assert inspect.isfunction(fix_Structure._structure_parser)


def test_has_datastr(fix_Structure):
    """Test has attribute."""
    assert isinstance(fix_Structure._datastr, str)


def test_filters_empty(fix_Structure):
    """Test self._filters is empty list."""
    assert not fix_Structure._filters


def test_kwargs_is_dict(fix_Structure):
    """Test **kwargs were stored to kwargs."""
    assert isinstance(fix_Structure.kwargs, dict)


def test_not_data_array(fix_Structure):
    """Test AttributeError."""
    with pytest.raises(EXCPTS.NotBuiltError):
        fix_Structure.data_array
# end


# # Test conditions after Strucure.build
def test_Structure_build_1(fix_Structure_build):
    """Build method."""
    with pytest.raises(AttributeError):
        fix_Structure_build._datastr


def test_Structure_build_2(fix_Structure_build):
    """Test data array is numpy array."""
    assert isinstance(fix_Structure_build.data_array, np.ndarray)
# end


# # Test parsing
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
    s = libstructure.Structure(cif_example, label_priority='label')
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
    assert len(s.chain_set) == len(expected)
    fasta = s.fasta
    assert fasta == expected
    assert next(iter(fasta.values())) == 'AYIAKQRQISFVKSHF'


def test_Structure_chain_set_CIF_auth():
    """Test chain reading from auth label."""
    s = libstructure.Structure(cif_example, label_priority='auth')
    s.build()
    assert len(s.chain_set) == 2


def test_Structure_build_3(fix_Structure_build):
    """Test shape of data_array."""
    expected = (545, len(libpdb.atom_slicers))
    assert fix_Structure_build.data_array.shape == expected


@pytest.mark.parametrize(
    'label_priority,expected',
    [
        ('label', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}),
        ('auth', {'A', 'B'}),
        ]
    )
def test_Structure_chainset_priority_label(label_priority, expected):
    """
    Test chain set with priority=label.

    This test is slighly repeated when testing `fasta`.
    However I found it appropriate to isolate it here.
    """
    s = libstructure.Structure(cif_example, label_priority=label_priority)
    s.build()
    assert s.chain_set == expected
# end


# Test filters
def test_false_num_of_filters(filtered):
    """Test number of filters."""
    assert len(filtered[0].filters) == 1


def test_false_len_filtered_atoms(filtered):
    """Test len of filtered atoms."""
    assert len(list(filtered[0].filtered_atoms)) == filtered[1]


def test_Structure_filter_chain_B(fix_Structure_build):
    """Test Structure filter chain B for CIF in `auth`."""
    # vim regex applied to count number of lines
    # on the pdb file
    # :%s/\(^ATOM\|^HETATM\).\{} B \s.//n
    # on the cif_file
    # %s/\(^ATOM\|^HETATM\).*\sB\s\D//n
    fix_Structure_build.add_filter_chain('B')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 275


def test_Structure_filter_B_cif_label():
    """Test Structure filter chain B for CIF in `label`."""
    # %s/\(^ATOM\|^HETATM\).*\sB\s\d//n
    s = libstructure.Structure(cif_example, label_priority='label')
    s.build()
    s.add_filter_chain('B')
    result = list(s.filtered_atoms)
    assert len(result) == 142


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


def test_Structure_filter_backbone_1(fix_Structure_build):
    """Test filter backbone atoms."""
    # CIF VIM REGEX: :%s/^ATOM.\{10}\(N\|CA\|C\|O\)\s.//n
    # PDB VIM REGEX: :%s/^ATOM.\{9}\(N\|CA\|C\|O\)\s.//n
    assert libstructure.col_name == 2
    fix_Structure_build.add_filter_backbone()
    fix_Structure_build.add_filter_record_name('ATOM')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 132


def test_Structure_filter_backbone_2(fix_Structure_build):
    """Test filter backbone atoms."""
    # CIF VIM REGEX: :%s/^ATOM.\{10}\(N\|CA\|C\|O\)\s.//n
    # PDB VIM REGEX: :%s/^ATOM.\{9}\(N\|CA\|C\|O\)\s.//n
    assert libstructure.col_name == 2
    fix_Structure_build.add_filter_backbone(minimal=True)
    fix_Structure_build.add_filter_record_name('ATOM')
    result = list(fix_Structure_build.filtered_atoms)
    assert len(result) == 132 // 4 * 3


def test_Structure_clear_filters(fix_Structure_build):
    """Test clear_filters method empties filters."""
    fix_Structure_build.add_filter(lambda x: True)
    fix_Structure_build.add_filter(lambda x: False)
    fix_Structure_build.add_filter(lambda x: None)
    assert len(fix_Structure_build.filters) == 3
    fix_Structure_build.clear_filters()
    assert len(fix_Structure_build.filters) == 0


def test_Structure_pop_last_filter(fix_Structure_build):
    """Test pop last filter."""
    def f1(): return True  # noqa: E704
    def f2(): return False  # noqa: E704
    def f3(): return 'a'  # noqa: E704
    fix_Structure_build.add_filter(f1)
    fix_Structure_build.add_filter(f2)
    fix_Structure_build.add_filter(f3)
    assert fix_Structure_build.filters[-1] is f3
    fix_Structure_build.pop_last_filter()
    assert fix_Structure_build.filters[-1] is f2
    fix_Structure_build.pop_last_filter()
    assert fix_Structure_build.filters[-1] is f1
    fix_Structure_build.pop_last_filter()
    assert len(fix_Structure_build.filters) == 0
# end test filters


# # Test methods
def test_consecutive_residues():
    """Test provides proper consecutive residues."""
    s = libstructure.Structure(pdb_res_gap)
    s.build()
    cr = s.consecutive_residues
    assert isinstance(cr, GeneratorType)
    result = [i for i in s.consecutive_residues]
    expected = [list(range(4, 11)), list(range(12, 20))]
    assert result == expected


def test_filtered_residues(fix_Structure_build):
    """Test Structure can list residues in filtered atoms."""
    fix_Structure_build.add_filter_chain('A')
    fix_Structure_build.add_filter(
        lambda x: int(x[libstructure.col_resSeq]) <= 10)
    result = fix_Structure_build.filtered_residues
    expected = list(range(4, 11))
    assert result == expected
# end


# # Test writing
def test_getPDB_with_filters(fix_Structure_build):
    """
    Test getPDB with filters.

    getPDB it self is tested in the following functions.
    """
    def line_filter(lines):
        for _ in lines:
            yield 80 * ' '
    lines = fix_Structure_build.get_PDB(pdb_filters=[line_filter])
    assert ''.join(lines).strip() == ''


def test_Structure_write(fix_Structure_build):
    """Test save functionality."""
    fix_Structure_build.add_filter_record_name('ATOM')
    fix_Structure_build.add_filter(
        lambda x: int(x[libstructure.col_resSeq]) < 11
        )
    fix_Structure_build.add_filter_chain('A')

    with TmpFile('pdb_testing.pdb') as fout:
        fix_Structure_build.write_PDB(fout)
        result = fout.read_text()
        expected = pdb_saved.read_text()
        assert result == expected


@pytest.mark.parametrize(
    'pdbmodels_in',
    [
        pdb_models,
        pdb_models2,
        ]
    )
def test_Structure_write_pdb_models(pdbmodels_in):
    """Test Structure parsig MODELS."""
    s = libstructure.Structure(pdbmodels_in)
    s.build()
    with TmpFile('pdb_testing.pdb') as fout:
        s.write_PDB(fout)
        result = fout.read_text()
        expected = pdb_models_output.read_text()
        assert result == expected


def test_Structure_write_empty_filter(fix_Structure_build):
    """Test save functionality."""
    fix_Structure_build.add_filter(lambda x: False)

    fout = Path('pdb_testing.pdb')
    with pytest.raises(EXCPTS.EmptyFilterError):
        fix_Structure_build.write_PDB(fout)


def test_function_writePDB(fix_Structure_build):
    """Test writePDB function sending manually from Structure."""
    fix_Structure_build.add_filter(lambda x: False)
    lines = libstructure.structure_to_pdb(fix_Structure_build.filtered_atoms)
    with pytest.warns(UserWarning):
        libstructure.write_PDB(lines, 'dummy_file.pdb')


def test_save_chains():
    """Test parse PDB in separate chains for IDPConfGen download pipeline."""
    execute = libstructure.save_structure_by_chains(
        pdb_example,
        'EXPL',
        chains=['A', 'B'],
        record_name=('ATOM',),
        )

    expected_fouts = ['EXPL_A.pdb', 'EXPL_B.pdb']
    expected_lines = [EXPL_A, EXPL_B]

    for i, (fout, lines) in enumerate(execute):
        assert fout == expected_fouts[i]
        explines = expected_lines[i].read_text().rstrip('\n')
        assert len(lines) == len(explines)
        # here, I have not been able to execute:
        # assert lines == explines
        # I can do it in the Jupyter notebook, but here CPU goes overwhelmed


def test_save_chain_ignored():
    """Test parse PDB in separate chains for IDPConfGen download pipeline."""
    execute = libstructure.save_structure_by_chains(
        pdb_example,
        '1MH1',
        chains=['A'],
        record_name=('ATOM',),
        )

    result = [i for i in execute]
    assert len(result) == 0


def test_save_chain_does_not_exist():
    """Test parse PDB in separate chains for IDPConfGen download pipeline."""
    execute = libstructure.save_structure_by_chains(
        pdb_example,
        'EXPL',
        chains=['Z'],
        record_name=('ATOM',),
        )

    result = [i for i in execute]
    assert len(result) == 0
# end


def test_generate_labels_1(fix_Structure_build):
    """Test generate labels function."""
    s = fix_Structure_build
    label_cols = [libstructure.col_resSeq, libstructure.col_resName, libstructure.col_name]
    labels = s.data_array[:3, label_cols]
    results = list(libstructure.concatenate_residue_labels(labels))
    expected = [
        '4ALAN',
        '4ALACA',
        '4ALAC',
        ]

    assert results == expected


def test_generate_labels_2(fix_Structure_build):
    """Test generate labels function."""
    s = fix_Structure_build
    label_cols = [libstructure.col_resSeq, libstructure.col_resName, libstructure.col_name]
    labels = s.data_array[:3, label_cols]
    labels1 = s.data_array[-3:, label_cols]
    results = libstructure.generate_residue_labels(
        labels,
        labels1,
        delimiter=',',
        )

    expected = [
        '4ALAN,431HOHO     ',
        '4ALACA,432HOHO    ',
        '4ALAC,433HOHO     ',
        ]

    assert results == expected

