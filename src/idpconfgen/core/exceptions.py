"""IDP Conf Gen Exceptions."""
from idpconfgen import contactus as CONTACTUS


class IDPConfGenException(Exception):
    """Base class exception for IDPConfGen."""

    errmsg = 'An error as occurred: {}. ' + CONTACTUS.contact_message

    def __init__(self, *args):
        self.args = args or ['']

    def __str__(self):
        """Represent exception as string."""
        return self.errmsg.format(*self.args)

    def __repr__(self):
        """Represent exception."""
        return f'{self.__class__.__name__}: {self}'


class PDBIDFactoryError(IDPConfGenException):
    """General PDBIDFactory Exception."""

    pass


class DownloadFailedError(IDPConfGenException):
    """Raise when download fails."""

    pass


class EmptyFilterError(IDPConfGenException):
    """Raise when PDB data filtering returns an empty selection."""
    errmsg = 'Filter returns empty selection.'


class DSSPParserError(IDPConfGenException):
    errmsg = 'Error while parsing {}'


class DSSPSecStructError(IDPConfGenException):
    errmsg = 'Values differ from possible DSSP secondary structure keys.'


class DSSPInputMissing(IDPConfGenException):
    errmsg = 'One of the two required positional arguments is missing.'