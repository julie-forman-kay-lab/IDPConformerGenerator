"""Test custom Exceptions."""

from idpconfgen import contactus as CONTACTUS
from idpconfgen.core import exceptions as EXCPTNS


class Test_General:
    """Test general Exception."""

    e = EXCPTNS.IDPConfGenException()

    def test_idpconfgen_general(self):
        """Test without message."""
        assert str(self.e) == \
            'An error as occurred: . ' + CONTACTUS.contact_message

    def test_repr(self):
        """Test repr."""
        assert repr(self.e) == (
            'IDPConfGenException: '
            'An error as occurred: . '
            f'{CONTACTUS.contact_message}'
            )