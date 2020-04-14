"""Test libcif."""
import pytest

from idpconfgen.core import exceptions as EXCPTS
from idpconfgen.libs.libcif import (
    CIFParser,
    is_cif,
    find_cif_atom_site_headers,
    populate_cif_dictionary,
    parse_cif_line,
    )

from .tcommons import (
    cif_example,
    cif_example_auth,
    cif_example_noasymid,
    )


@pytest.mark.parametrize(
    'in1,expected',
    [
        (
            "ATOM   1    N N     . ALA A 1 4   ? 11.751 37.846 29.016  1.00 44.65 ? 4   ALA A N     1 \n",
            ['ATOM', '1', 'N', 'N', '.', 'ALA', 'A', '1', '4', '?', '11.751', '37.846', '29.016', '1.00', '44.65', '?', '4', 'ALA', 'A', 'N', '1'],
            )
        ]
    )
def test_cif_line_regex(in1, expected):
    """Test cif line parse according to regex."""
    result = parse_cif_line(in1)
    assert result == expected


def test_is_cif():
    """Test file is cif."""
    assert is_cif(cif_example.read_text())


@pytest.fixture(params=[cif_example])
def fix_find_cif_atom_site_headers(request):
    """Execute find CIF atom site headers."""
    lines = cif_example.read_text().split('\n')
    cif_pdb = {}
    index = find_cif_atom_site_headers(lines, cif_pdb)
    return index, cif_pdb, lines


def test_find_cif_atom_site_headers_index(fix_find_cif_atom_site_headers):
    """Text index matches first structure line after atom_site headers."""
    index = fix_find_cif_atom_site_headers[0]
    assert index == 1192


def test_find_cif_atom_site_headers_len(fix_find_cif_atom_site_headers):
    """Test the number of headings found."""
    cif_pdb = fix_find_cif_atom_site_headers[1]
    assert len(cif_pdb) == 21


def test_find_cif_atom_site_headers_keys(fix_find_cif_atom_site_headers):
    """Test all headings are `_atom_site.`"""
    cif_pdb = fix_find_cif_atom_site_headers[1]
    assert all(k.startswith('_atom_site.') for k in cif_pdb.keys())


def test_fine_cif_atom_site_headers_error():
    """Test missing atom_site headers error."""
    with pytest.raises(EXCPTS.CIFFileError):
        find_cif_atom_site_headers([' ', 'somestring'], {})


def test_populate_cif_dictionary(fix_find_cif_atom_site_headers):
    """Test populate cif dictionary."""
    index = fix_find_cif_atom_site_headers[0]
    cif_pdb = fix_find_cif_atom_site_headers[1]
    lines = fix_find_cif_atom_site_headers[2]

    number_of_atoms = populate_cif_dictionary(lines, index, cif_pdb)
    assert number_of_atoms == 5385
    assert cif_pdb['_atom_site.label_atom_id'][:5] == \
        ['N', 'CA', 'C', 'O', 'CB']
    assert cif_pdb['_atom_site.Cartn_x'][-2] == '28.365'


@pytest.fixture(params=[
    cif_example,
    cif_example_auth,
    ])
def cif_parsed(request):
    return CIFParser(request.param.read_text())


def test_CIFParse_len(cif_parsed):
    """Test CIFParse len."""
    assert len(cif_parsed) == 5385


def test_CIFParse_line(cif_parsed):
    """Test CIFParse len."""
    cif_parsed.line = 5
    assert cif_parsed.line == 5


def test_CIFParse_line_error(cif_parsed):
    with pytest.raises(AssertionError):
        cif_parsed.line = -5

def test_CIFParse_get_values(cif_parsed):
    cif_parsed.line = 0
    expected = [
        'ATOM',
        '1',
        'N',
        ' ',
        'ALA',
        'A',
        '4',
        ' ',
        '11.751',
        '37.846',
        '29.016',
        '1.00',
        '44.65',
        ' ',
        'N',
        ' ']
    assert cif_parsed.get_line_elements_for_PDB() == expected


def test_CIFParse_get_values_line_2(cif_parsed):
    expected = [
        'ATOM',
        '3',
        'C',
        ' ',
        'ALA',
        'A',
        '4',
        ' ',
        '13.740',
        '38.628',
        '27.754',
        '1.00',
        '24.74',
        ' ',
        'C',
        ' ']
    assert cif_parsed.get_line_elements_for_PDB(line=2) == expected


def test_CIFParse_get_values_no_chaindi():
    """
    Also test lack of ins code.
    """
    cif = CIFParser(cif_example_noasymid.read_text())
    expected = [
        'ATOM',
        '3',
        'C',
        ' ',
        'ALA',
        'A',
        '4',
        ' ',
        '13.740',
        '38.628',
        '27.754',
        '1.00',
        '24.74',
        ' ',
        'C',
        ' ']
    assert cif.get_line_elements_for_PDB(line=2) == expected
