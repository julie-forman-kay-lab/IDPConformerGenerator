"""Test libcli."""

from idpconfgen.libs import libcli


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
            '\n'
            '    $ use the program this way\n'
            '    $ you can also use it like this.\n'
            )
        assert expected == self.usage
