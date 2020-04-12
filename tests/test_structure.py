"""Test Structure class."""
import pytest
import hypothesis.strategies as st
from hypothesis import given


from idpconfgen.libs.libpdb import (
    gen_pdb_data_array,
    is_cif,
    is_pdb,
    parse_pdb_to_array,
    parse_cif_to_array,
    PDBParams,
    CIFParser,
    Structure
    )

from . tcommons import (
    pdb_example,
    cif_example,
    )


def test_PDBParams():
    """
    """
    with pytest.raises(NotImplementedError):
        PDBParams.some_attr = None


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
    result = gen_pdb_data_array(number)
    assert result.shape == (number, len(PDBParams.atom_slicers))
    assert result.dype == '<U8'


def test_parse_pdb_to_array():
    """
    """
    data_array = parse_pdb_to_array(pdb_example.read_text())
    assert data_array.shape == (3000, len(PDBParams.atom_slicers))
    assert data_array.dtype == np.dtype('<U8')


# TODO
def test_parse_cif_to_array():
    """
    """
    data_array = parse_cif_to_array(cif_example.read_text())
    assert data_array.shape == (3000, len(PDBParams.atom_slicers))
    assert data_array.dtype == np.dtype('<U8')


@pytest.fixture()
def CIFParser_fixture(request):
    return CIFParser(cif_example.read_text())

@pytest.mark.parametrize(
    'line,expected',
    [
        # TODO
        ]
    )
def test_CIFParser_init(CIFParser_fixture, line, expected):
    """
    """
    cif = CIFParser_fixture
    assert cif.get_line_elements_for_PDB(line) == expected


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


@pytest.fixture
def Structure_built(request):
    s = Structure(pdb_example.read_text())
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
