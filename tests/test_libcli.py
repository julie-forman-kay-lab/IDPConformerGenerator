"""Test libcli."""
import argparse
import os
import sys

import pytest

from idpconfgen import Path
from idpconfgen.libs import libcli


def test_load_args():
    """Test load args."""
    ap = argparse.ArgumentParser()
    ap.add_argument('dummy')
    prev = sys.argv
    sys.argv = ['blab', 'this_is_dummy']
    result = libcli.load_args(ap)
    expected = ap.parse_args(['this_is_dummy'])
    assert expected == result
    sys.argv = prev


def test_maincli():
    """Test maincli."""
    def main(a, b=1):
        assert a == 'namea'
        assert b == 9
        return
    ap = argparse.ArgumentParser()
    ap.add_argument('a')
    ap.add_argument('-b', type=int)
    prev = sys.argv
    args = ['namea', '-b', '9']
    sys.argv = ['mycommand'] + args
    libcli.maincli(ap, main)
    sys.argv = prev


@pytest.mark.parametrize(
    'dest,expected',
    [
        ('myfolder', ['myfolder']),
        ('mytar.tar', 'mytar.tar'),
        ]
    )
def test_folderortar(dest, expected):
    """Test folder or tar."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        'dest',
        action=libcli.FolderOrTar,
        nargs='+',
        )
    result = ap.parse_args([dest])
    assert result.dest == expected


@pytest.mark.parametrize(
    'args,expected',
    [
        (['a', 'b', 'c'], ('a', 'b', 'c')),
        (['abc'], ('abc',)),
        (['abc', '1', 'z'], ('abc', '1', 'z')),
        (['ATOM'], ('ATOM',)),
        (['ATOM', 'HETATM'], ('ATOM', 'HETATM')),
        ]
    )
def test_args2tuple(args, expected):
    """Convert args to tuple."""
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', action=libcli.ArgsToTuple, nargs='+')
    result = ap.parse_args(['-t'] + args)
    print(result)
    assert result.t == expected


@pytest.mark.parametrize(
    'ext, args',
    [
        ['.json', 'mydict.json'],
        ['.tar', 'myfile.tar'],
        ]
    )
def test_checkext(ext, args):
    """Check Extensions."""
    ap = argparse.ArgumentParser()
    ap.add_argument('-e', action=libcli.CheckExt(ext))
    result = ap.parse_args(['-e'] + [args])
    assert result.e == args


def test_SeqOrFasta(in1, expected):


def test_checkext_error():
    """Check Ext with error."""
    ap = argparse.ArgumentParser()
    ap.add_argument('-e', action=libcli.CheckExt('.json'))
    with pytest.raises(SystemExit):
        ap.parse_args(['-e'] + ['myfile.txt'])


# capsys, crazy pytest undeclared variable! :-o
# https://docs.pytest.org/en/stable/reference.html?highlight=capsys#capsys
def test_customparse(capsys):
    """Test custom parsers."""
    assert issubclass(libcli.CustomParser, argparse.ArgumentParser)
    ap = libcli.CustomParser()
    with pytest.raises(SystemExit) as err:
        ap.error('message')
    assert err.value.args[0] == 2
    captured = capsys.readouterr()
    assert captured.err == '\nerror: message\n'


def test_add_subparse():
    """Test add subparser."""
    ap = argparse.ArgumentParser(description='description')
    nspace = argparse.Namespace()
    nspace._prog = 'somecli'
    nspace._name = 'mycli'
    nspace._help = 'cli help message'
    nspace._usage = 'myusage'
    nspace.main = lambda x: True
    nspace.ap = ap

    ap2 = argparse.ArgumentParser()
    subparser = ap2.add_subparsers()
    libcli.add_subparser(subparser, nspace)


def test_add_version():
    """Test add version."""
    ap = argparse.ArgumentParser()
    libcli.add_version(ap)
    with pytest.raises(SystemExit) as err:
        ap.parse_args(['-v'])
    assert err.value.args[0] == 0


@pytest.mark.parametrize(
    'func, parsing, expected',
    [
        (libcli.add_argument_chunks, ['-c', '5000'], ['chunks', 5000]),
        (libcli.add_argument_destination_folder, ['-d', 'myfolder'], ['destination', Path('myfolder')]),  # noqa: E501
        (libcli.add_argument_minimum, ['-m', '45'], ['minimum', 45]),
        (libcli.add_argument_ncores, ['-n', '7'], ['ncores', 7]),
        (libcli.add_argument_ncores, ['-n'], ['ncores', os.cpu_count() - 1]),
        (libcli.add_argument_output, [], ['output', None]),
        (libcli.add_argument_output, ['-o', 'out.log'], ['output', 'out.log']),
        (libcli.add_argument_output, ['--output', 'out.log'], ['output', 'out.log']),
        (libcli.add_argument_pdb_files, ['myfolder'], ['pdb_files', ['myfolder']]),  # noqa: E501
        (libcli.add_argument_pdb_files, ['my.tar'], ['pdb_files', 'my.tar']),
        (libcli.add_argument_pdbids, ['12AS_A'], ['pdbids', ['12AS_A']]),
        (libcli.add_argument_reduced, ['-rd'], ['reduced', True]),
        (libcli.add_argument_reduced, [], ['reduced', False]),
        (libcli.add_argument_replace, ['--replace'], ['replace', True]),
        (libcli.add_argument_replace, [], ['replace', False]),
        (libcli.add_argument_record, ['-rn', 'ATOM'], ['record_name', ('ATOM',)]),  # noqa: E501
        (libcli.add_argument_source, ['-sc', 'source.txt'], ['source', Path('source.txt')]),  # noqa: E501
        (libcli.add_argument_update, ['-u'], ['update', True]),
        (libcli.add_argument_update, [], ['update', False]),
        (libcli.add_argument_cmd, ['mycmd'], ['cmd', 'mycmd']),
        ]
    )
def test_add_arguments(func, parsing, expected):
    """Test add arguments."""
    ap = argparse.ArgumentParser()
    func(ap)
    nspace = ap.parse_args(parsing)
    assert getattr(nspace, expected[0]) == expected[1]


class TestParseDoc:
    """Test parsing docstring from cli scripts."""

    docstring = """
Program Name

DESCRIPTION:
    
    This is the description of the program.
    In two lines.

USAGE:

    $ use the program this way
    $ you can also use it like this.
"""
    
    prog, des, usage = libcli.parse_doc_params(docstring)

    def test_prog(self):
        """Test prog string description."""
        expected_prog = "Program Name"
        assert expected_prog == self.prog

    def test_description(self):
        """Test description string description."""
        expected = (
            'DESCRIPTION:\n'
            '    \n'
            '    This is the description of the program.\n'
            '    In two lines.\n'
            )
        print(self.des)
        assert expected == self.des

    def test_usage(self):
        """Test usage string description."""
        expected = (
            '\n\n'
            '    $ use the program this way\n'
            '    $ you can also use it like this.\n'
            )
        assert expected == self.usage
