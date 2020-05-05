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
    cif_nohash,
    cif_example_headers,
    cif_EOF,
    )


@pytest.mark.parametrize(
    'in1,expected',
    [
        (
            "ATOM   1    N N     . ALA A 1 4   ? 11.751 37.846 29.016  1.00 44.65 ? 4   ALA A N     1 \n",
            ['ATOM', '1', 'N', 'N', '.', 'ALA', 'A', '1', '4', '?', '11.751', '37.846', '29.016', '1.00', '44.65', '?', '4', 'ALA', 'A', 'N', '1'],
            ),
        (
            "HETATM 5119 N N     . ASN C 2 .   ? 29.722 14.604 9.802   1.00 12.15 ? 331 ASN A N     1 ",
            ['HETATM', '5119', 'N', 'N', '.', 'ASN', 'C', '2', '.', '?', '29.722', '14.604', '9.802', '1.00', '12.15', '?', '331', 'ASN', 'A', 'N', '1'],
            ),
        (
            'HETATM 5132 O "O5\'" . AMP D 3 .   ? 23.706 19.236 9.885   1.00 10.43 ? 332 AMP A "O5\'" 1 ',
            ['HETATM', '5132', 'O', "O5'", '.', 'AMP', 'D', '3', '.', '?', '23.706', '19.236', '9.885', '1.00', '10.43', '?', '332', 'AMP', 'A', "O5'", '1'],
            ),
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
    assert index == 58


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
    assert number_of_atoms == 545
    assert cif_pdb['_atom_site.label_atom_id'][:5] == \
        ['N', 'CA', 'C', 'O', 'CB']
    assert cif_pdb['_atom_site.Cartn_x'][-2] == '28.365'


@pytest.mark.parametrize(
    'args',
    [
        (cif_EOF.read_text().split('\n'), 22, {i:[] for i in cif_example_headers}),
        (cif_nohash.read_text().split('\n'), 22, {i:[] for i in cif_example_headers}),
        ([''], 0, {'_atom_site.id': []}),
        ]
    )
def test_populate_cif_dictionary_errors(args):
    """Test CIF without hash raise error."""
    with pytest.raises(EXCPTS.CIFFileError) as err:
        populate_cif_dictionary(*args)


@pytest.fixture(params=[
    cif_example,
    cif_example_auth,
    ])
def cif_parsed(request):
    return CIFParser(request.param.read_text())


def test_CIFParse_len(cif_parsed):
    """Test CIFParse len."""
    assert len(cif_parsed) == 545

def test_CIFParse_line_default(cif_parsed):
    """Test default attribute line."""
    assert cif_parsed.line == 0

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


def test_CIFParse_get_values_line_3(cif_parsed):
    expected = [
        'HETATM',
        '5385',
        'O',
        ' ',
        'HOH',
        'H',
        ' ',
        ' ',
        '45.659',
        '28.060',
        '40.225',
        '1.00',
        '3.27',
        ' ',
        'O',
        ' ']
    assert cif_parsed.get_line_elements_for_PDB(line=544) == expected



def test_CIFParse_get_values_no_chainid():
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
