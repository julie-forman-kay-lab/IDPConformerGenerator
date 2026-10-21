"""Test support for DSSP v2, v3 and v4."""
import shlex
import shutil

import pytest

from idpconfgen import Path
from idpconfgen.cli_dssppii import assign_ppii, dssp_ppii_assignment
from idpconfgen.cli_sscalc import dssppi_helper
from idpconfgen.core.exceptions import DSSPExecutionError, DSSPParserError
from idpconfgen.libs import libdssp
from idpconfgen.libs.libparse import split_pdb_by_dssp

from . import tcommons


PDB_1BTE = Path(tcommons.data_folder, '1BTE_A.pdb')

# Output of `mkdssp` for 1BTE_A, by DSSP version
DSSP_1BTE = {
    '3': Path(tcommons.data_folder, '1BTE_A.dssp'),  # DSSP 3.1.2
    '4': Path(tcommons.data_folder, '1BTE_A_dssp4.dssp'),  # DSSP 4.6.1
    }

# Scripts that emulate the command line of `mkdssp` v3 and v4. The v4 script
# also fails, as DSSP 4.2.2 does, on PDB files without a HEADER line at the
# top (DSSP 4.2.2 reads them as mmCIF) or without a newline at the end
# (DSSP 4.2.2 never finishes reading them).
FAKE_MKDSSP = {
    '3': (
        '#!/bin/sh\n'
        'if [ "$1" = "--version" ]; then\n'
        '    echo "mkdssp version 3.1.4"\n'
        '    echo\n'
        '    echo "If you use DSSP, please cite: "\n'
        'elif [ "$1" = "-i" ]; then\n'
        '    cat {output}\n'
        'else\n'
        '    echo "unrecognised option" >&2\n'
        '    exit 1\n'
        'fi\n'
        ),
    '4': (
        '#!/bin/sh\n'
        'if [ "$1" = "--version" ]; then\n'
        '    echo "mkdssp version 4.2.2"\n'
        'elif [ "$1" = "--output-format" ] && [ "$2" = "dssp" ]; then\n'
        '    if [ -n "$(tail -c 1 "$3")" ]; then\n'
        '        echo "no newline at the end of the file" >&2\n'
        '        exit 1\n'
        '    fi\n'
        '    read -r first_line < "$3"\n'
        '    case "$first_line" in\n'
        '        HEADER*) cat {output} ;;\n'
        '        *) echo "not an mmCIF file" >&2; exit 1 ;;\n'
        '    esac\n'
        'else\n'
        '    echo "unknown option" >&2\n'
        '    exit 1\n'
        'fi\n'
        ),
    }


@pytest.fixture
def fake_mkdssp(tmp_path):
    """
    Return a function that makes a fake `mkdssp` executable.

    The fake `mkdssp` accepts the command line of DSSP v3 or v4 and
    prints the output of that DSSP version for 1BTE_A.
    """
    def make(version):
        exe = Path(tmp_path, f'mkdssp{version}')
        output = shlex.quote(str(DSSP_1BTE[version]))
        exe.write_text(FAKE_MKDSSP[version].format(output=output))
        exe.chmod(0o755)
        return str(exe)

    return make


def _sscalc_records(dssp_lines):
    """Split 1BTE_A with DSSP data as `idpconfgen sscalc` does."""
    dssp_text = '\n'.join(dssp_lines).encode()
    return [
        (name, data)
        for name, data, _ in split_pdb_by_dssp(PDB_1BTE, dssp_text)
        ]


@pytest.fixture
def dssp3_records():
    """Return the `sscalc` data for 1BTE_A from DSSP v3 with DSSP-PPII."""
    lines = assign_ppii(DSSP_1BTE['3'].read_text().splitlines())
    return _sscalc_records(lines)


@pytest.mark.parametrize(
    'text, expected',
    [
        ('mkdssp version 3.1.4\n\nIf you use DSSP, please cite:\n', (3, 1, 4)),
        ('mkdssp version 4.6.1\n', (4, 6, 1)),
        ('mkdssp version 2.2.1\n', (2, 2, 1)),
        ('mkdssp version 4.2.2.1\n', (4, 2, 2, 1)),
        ('/usr/bin/mkdssp version 4.0.0\n', (4, 0, 0)),
        ('Usage: mkdssp [options] input-file [output-file]\n', None),
        ('', None),
        ]
    )
def test_parse_dssp_version(text, expected):
    """Test reading the DSSP version."""
    assert libdssp.parse_dssp_version(text) == expected


@pytest.mark.parametrize(
    'version, expected, as_str',
    [
        ((2, 2, 1), False, '2.2.1'),
        ((3, 1, 4), False, '3.1.4'),
        ((4, 0, 4), True, '4.0.4'),
        ((4, 6, 1), True, '4.6.1'),
        (None, False, 'unknown'),
        ]
    )
def test_is_dssp4(version, expected, as_str):
    """Test DSSP v4 is identified from its version."""
    assert libdssp.is_dssp4(version) is expected
    assert libdssp.dssp_version_str(version) == as_str


@pytest.mark.parametrize(
    'version, expected',
    [
        ((3, 1, 4), ['mkdssp', '-i', 'x.pdb']),
        ((4, 6, 1), ['mkdssp', '--output-format', 'dssp', 'x.pdb']),
        (None, ['mkdssp', '-i', 'x.pdb']),
        ]
    )
def test_make_dssp_args(version, expected):
    """Test DSSP is called according to its version."""
    assert libdssp.make_dssp_args('x.pdb', 'mkdssp', version) == expected


def test_split_dssp_cmd():
    """Test commands with arguments are split."""
    result = libdssp.split_dssp_cmd('singularity exec "my dssp.sif" mkdssp')
    assert result == ['singularity', 'exec', 'my dssp.sif', 'mkdssp']


def test_split_dssp_cmd_path_with_spaces(tmp_path):
    """Test paths to an executable are not split at spaces."""
    exe = Path(tmp_path, 'my dssp', 'mkdssp')
    exe.parent.mkdir()
    exe.write_text('#!/bin/sh\n')
    exe.chmod(0o755)
    assert libdssp.split_dssp_cmd(str(exe)) == [str(exe)]


def test_split_dssp_cmd_empty():
    """Test empty DSSP commands raise an error."""
    with pytest.raises(DSSPExecutionError):
        libdssp.split_dssp_cmd('')


@pytest.mark.parametrize(
    'version, expected',
    [
        ('3', (3, 1, 4)),
        ('4', (4, 2, 2)),
        ]
    )
def test_get_dssp_version(fake_mkdssp, version, expected):
    """Test reading the version of the DSSP executable."""
    assert libdssp.get_dssp_version(fake_mkdssp(version)) == expected


def test_get_dssp_version_not_found():
    """Test missing DSSP executables raise an error."""
    with pytest.raises(DSSPExecutionError):
        libdssp.get_dssp_version('idpconfgen-no-such-mkdssp')


@pytest.mark.parametrize('version', ['3', '4'])
def test_run_dssp(fake_mkdssp, version):
    """Test DSSP v3 and v4 are called with arguments they accept."""
    result = libdssp.run_dssp(PDB_1BTE, fake_mkdssp(version))
    assert result == DSSP_1BTE[version].read_bytes()


def test_run_dssp_no_final_newline(fake_mkdssp, tmp_path):
    """Test DSSP v4 reads PDB files whose last line has no newline."""
    pdb = Path(tmp_path, '1BTE_A.pdb')
    pdb.write_bytes(PDB_1BTE.read_bytes().rstrip(b'\n'))
    result = libdssp.run_dssp(pdb, fake_mkdssp('4'))
    assert result == DSSP_1BTE['4'].read_bytes()


def test_run_dssp_error(fake_mkdssp):
    """Test DSSP errors are reported, here calling DSSP v4 as DSSP v3."""
    with pytest.raises(DSSPExecutionError) as error:
        libdssp.run_dssp(PDB_1BTE, fake_mkdssp('4'), dssp_version=(3, 1, 4))
    assert 'exited with code 1' in str(error.value)
    assert 'unknown option' in str(error.value)


def test_run_dssp_killed(tmp_path):
    """Test DSSP killed by a signal is reported with the signal name."""
    exe = Path(tmp_path, 'mkdssp')
    exe.write_text('#!/bin/sh\nkill -9 $$\n')
    exe.chmod(0o755)
    with pytest.raises(DSSPExecutionError) as error:
        libdssp.run_dssp(PDB_1BTE, str(exe), dssp_version=(4, 2, 2))
    assert 'was killed by SIGKILL' in str(error.value)


@pytest.mark.parametrize(
    'top, bottom',
    [
        (b'', b''),  # as IDPConformerGenerator writes PDB files
        (b'', b'\n'),
        (b'HEADER    TEST\n', b''),
        ]
    )
def test_pdb_for_dssp4(tmp_path, top, bottom):
    """Test DSSP v4 gets a copy with a HEADER line and a final newline."""
    atoms = PDB_1BTE.read_bytes().rstrip(b'\n')
    pdb = Path(tmp_path, 'test.pdb')
    pdb.write_bytes(top + atoms + bottom)

    with libdssp.pdb_for_dssp4(pdb) as dssp_input:
        tmp_pdb = Path(dssp_input)
        assert tmp_pdb != pdb
        data = tmp_pdb.read_bytes()

    assert data == (top or libdssp.PDB_HEADER_LINE) + atoms + b'\n'
    assert not tmp_pdb.exists()
    assert pdb.read_bytes() == top + atoms + bottom


def test_pdb_for_dssp4_no_copy(tmp_path):
    """Test PDB files with a HEADER line and a final newline are not copied."""
    pdb = Path(tmp_path, 'header.pdb')
    pdb.write_bytes(b'HEADER    TEST\n' + PDB_1BTE.read_bytes())
    with libdssp.pdb_for_dssp4(pdb) as dssp_input:
        assert dssp_input == pdb


def test_check_dssp_output():
    """Test classic DSSP output passes the check."""
    for dssp_file in DSSP_1BTE.values():
        libdssp.check_dssp_output(dssp_file.read_bytes())


@pytest.mark.parametrize(
    'output',
    [
        b'data_nohd\n# \n_entry.id   nohd\n# \n',  # mmCIF, from DSSP v4
        b'',
        ]
    )
def test_check_dssp_output_error(output):
    """Test output not in the classic DSSP format raises an error."""
    with pytest.raises(DSSPParserError):
        libdssp.check_dssp_output(output)


def test_remove_native_ppii():
    """Test the PPII helix of DSSP v4 in 1BTE_A is reset to loop."""
    lines = DSSP_1BTE['4'].read_text().splitlines()
    result = libdssp.remove_native_ppii(lines)
    assert len(result) == len(lines)

    changed = [i for i, line in enumerate(lines) if line != result[i]]
    assert len(changed) == 1  # Cys31 is the only residue with P

    old, new = lines[changed[0]], result[changed[0]]
    assert old[5:10].strip() == '31' and old[16] == 'P'
    assert new == old[:16] + ' ' + old[17:]


def test_remove_native_ppii_residue_table_only():
    """Test lines before the residue table are not changed."""
    lines = [
        'P' * 40,
        libdssp.DSSP_RESIDUE_HEADER,
        'x' * 16 + 'P' + 'x' * 10,
        ]
    result = libdssp.remove_native_ppii(lines)
    assert result == lines[:2] + ['x' * 16 + ' ' + 'x' * 10]


@pytest.mark.parametrize('version', ['3', '4'])
def test_dssp_ppii_assignment(fake_mkdssp, dssp3_records, version):
    """Test DSSP v3 and v4 give the same result with DSSP-PPII."""
    lines = dssp_ppii_assignment(str(PDB_1BTE), fake_mkdssp(version))
    records = _sscalc_records(lines)
    assert [name for name, _ in records] == ['1BTE_A_seg0', '1BTE_A_seg1']
    assert records == dssp3_records


def test_dssp_ppii_assignment_native_ppii(fake_mkdssp, dssp3_records):
    """Test `native_ppii` keeps the PPII helices assigned by DSSP v4."""
    lines = dssp_ppii_assignment(
        str(PDB_1BTE),
        fake_mkdssp('4'),
        native_ppii=True,
        )
    assert lines == DSSP_1BTE['4'].read_text().splitlines()

    # DSSP v4 assigns P to Cys31, DSSP-PPII does not
    seg0 = _sscalc_records(lines)[0][1]
    cys31 = seg0['resids'].split(',').index('31')
    assert seg0['dssp'][cys31] == 'P'
    assert dssp3_records[0][1]['dssp'][cys31] == ' '


@pytest.mark.parametrize('version', ['3', '4'])
def test_sscalc_helper(fake_mkdssp, dssp3_records, version):
    """Test `idpconfgen sscalc` gives the same data for DSSP v3 and v4."""
    exe = fake_mkdssp(version)
    results = dssppi_helper(
        PDB_1BTE,
        exe,
        dssp_version=libdssp.get_dssp_version(exe),
        minimum=2,
        reduced=False,
        )
    assert [(name, data) for name, data, _ in results] == dssp3_records


@pytest.mark.skipif(not shutil.which('mkdssp'), reason="requires mkdssp")
def test_dssp_ppii_assignment_mkdssp(dssp3_records):
    """Test the installed `mkdssp`, which can be DSSP v2, v3 or v4."""
    lines = dssp_ppii_assignment(str(PDB_1BTE), 'mkdssp')
    assert _sscalc_records(lines) == dssp3_records
