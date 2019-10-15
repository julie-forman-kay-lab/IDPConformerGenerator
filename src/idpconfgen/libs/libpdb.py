"""Contain  handlers of PDB information."""
import functools
from pathlib import Path
import re

from idpconfgen import log
from idpconfgen.core import exceptions as EXCPTS


class PDBParams:
    """
    PDB Format string slicing according to:
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    """
    line_name = slice(0, 6)
    serial = slice(6, 11)
    atom_name = slice(12, 16)
    altloc = slice(16, 17)
    resname = slice(17, 20)
    chainid = slice(21, 22)
    resseq = slice(22, 26)
    icode = slice(26, 27)
    xcoord = slice(30, 38)
    ycoord = slice(38, 46)
    zcoord = slice(46, 54)
    occupancy = slice(54, 60)
    tempfactor = slice(60, 66)
    element = slice(76, 78)
    charge = slice(78, 80)


class PDBIDFactory:
    """
    Parses input for PDBID instatiation.
    
    Returns a PDBID object.
    
    (okay... this might not be exactly a factory,
    but I wanted to have a object instantiation interface type
    instead of function call)
    """
    cull_onlypdbid_regex_single = re.compile(r'^[0-9a-zA-Z]{4}$')
    cull_onlypdbid_regex = re.compile(r'^[0-9a-zA-Z]{4}[\s\n]')
    cull_with_chains_regex_single = re.compile(r'^[0-9a-zA-Z]{5,}$')
    cull_with_chains_regex = re.compile(r'^[0-9a-zA-Z]{5,}[\s\n]')
    pdb_filename_regex = re.compile(r'[0-9a-zA-Z]{4}_[0-9a-zA-Z]{1,}\.pdb$')
    
    def __new__(cls, name):
       
        if isinstance(name, PDBID):
            return name
        
        pdb_filename_regex = {
            cls.cull_onlypdbid_regex_single: cls._parse_cull_onlypdbid,
            cls.cull_onlypdbid_regex: cls._parse_cull_onlypdbid,
            cls.cull_with_chains_regex_single: cls._parse_cull_with_chains,
            cls.cull_with_chains_regex: cls._parse_cull_with_chains,
            cls.pdb_filename_regex: cls._parse_filename_with_chain,
            }
        
        for regex, parser in pdb_filename_regex.items():
            if regex.search(str(name)):  # in case Path obj
                break
        else:
            emsg = f"PDB code format not valid: {name}. No regex matches."
            raise EXCPTS.PDBIDFactoryError(emsg)
        
        return PDBID(*parser(name))
    
    @staticmethod
    def _parse_cull_onlypdbid(pdbid):
        return pdbid[:4], None
    
    @staticmethod
    def _parse_cull_with_chains(pdbid):
        pdbinfo = pdbid.split()[0]
        return pdbinfo[:4], pdbinfo[4:]
    
    @staticmethod
    def _parse_filename_with_chain(pdbid):
        pdbid, chainid = Path(pdbid).stem.split('_')
        return pdbid, chainid


class PDBList:
    """
    List of PDBID objects.
    
    Parameters
    ----------
    pdb_names : obj:iterator
        An iterator containing the PDB names.
        PDB names can be in the form accepted by PDBIDFactory or
        PDBID objects.
    """
    def __new__(cls, pdb_names):
        
        try:
            if isinstance(pdb_names, cls):
                return pdb_names
            else:
                return super().__new__(cls)
        except IndexError:
            return super().__new__(cls)
        
    def __init__(self, pdb_names):
        valid_pdb_names = filter(
            lambda x: not str(x).startswith('#'),  # may receive Paths
            pdb_names,
            )
        self.set = set(PDBIDFactory(element) for element in valid_pdb_names)
    
    def __repr__(self):
        return '{}(\n    {})\n'.format(
            self.__class__.__name__,
            ',\n    '.join(repr(x) for x in self),
            )
    
    def __str__(self):
        return '{} with {} elements'.format(
            self.__class__.__name__,
            len(self),
            )
    
    def __eq__(self, other):
        try:
            return self.set == other.set
        except AttributeError:
            return self.set == other
    
    def __iter__(self):
        return iter(self.to_tuple())
    
    def __getitem__(self, index):
        return self.to_tuple()[index]
    
    def __len__(self):
        return len(self.set)
    
    def to_tuple(self):
        """
        Converts PDBList to sorted tuple.
        """
        return tuple(sorted(self.set))
    
    def difference(self, other):
        """
        Returns a PDBList that's the set difference between
        self and other.
        """
        return PDBList(tuple(self.set.difference(other.set)))
    
    def write(self, filename='PDBIDs.list'):
        """
        Writes to filename the PDBIDs in the PDBList.
        """
        with open(filename, 'w') as fh:
            fh.write('\n'.join(str(pdbid) for pdbid in self.to_tuple()))
        
        log.info(S(f'PDBIDs written to {filename}'))
        

@functools.total_ordering
class PDBID:
    """
    PDB object identifier.
    
    Identifies unique downloadable/stored units.
    
    In the current implmentation each unit is one PDB chain,
    which is identified by the PDBID and its chain identifier.
    
    Parameters
    ----------
    name : obj:`str`
        The PDBID, for example: 1ABC
    
    chain : obj:`str`
        The chain identifier.
        Defaults to None.
    
    Attributes
    ----------
    name:
        The four character PDB identifier.
    
    chain:
        The chain identifier.
    """
    def __init__(self, name, chain=None):
        
        self.name = name.upper()
        self.chain = chain
        
        # made manual to completely control order
        ids = {
            'chain': chain,
            }
        self.identifiers = {}
        for name, identifier in ids.items():
            if identifier:
                self.identifiers[name] = identifier
    
    def __repr__(self):
        
        iditems = self.identifiers.items()
        
        kwargs = ', '.join(f'{key}={val!r}' for key, val in iditems)
        
        if kwargs:
            kwargs = ', ' + kwargs
        
        return '{}(name={!r}{})'.format(
            self.__class__.__name__,
            self.name,
            kwargs,
            )
    
    def __lt__(self, other):
        return str(self) < str(other)
    
    def __hash__(self):
        return hash(str(self))
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    def __str__(self):
        name = f'{self.name}'
        ids = '_'.join(self.identifiers.values())
        
        if ids:
            return f'{name}_' + ids
        else:
            return name



class PDBDownloader:
    def __init__(
            self,
            pdb_list,
            destination,
            ncores=1,
            record_name=('ATOM',),
            **kwargs,
            ):
        """
        
        destination (folder)
        """
        
        self.pdb_list = pdb_list

        self.destination = Path(destination)
        self.destination.mkdir(parents=True, exist_ok=True)
        
        self.ncores = ncores
        
        self.record_name = record_name
        
        self.kwargs = kwargs
        
        self.prepare_pdb_list()
        
        return
    
    def prepare_pdb_list(self):
        
        self.pdbs_to_download = {}
        for pdbid in self.pdb_list:
            pdbentry = self.pdbs_to_download.setdefault(pdbid.name, [])
            pdbentry.append(pdbid.chain)
    
    @record_time()
    def run(self):
        log.info(T('starting raw PDB download'))
        log.info(S(
            f'{len(self.pdbs_to_download)} PDBs will be downloaded '
            f'and at least {len(self.pdb_list)} chains will be saved'
            ))
        
        results = ThreadPool(self.ncores).imap_unordered(
            self._download_single_pdb,
            self.pdbs_to_download.items(),
            )
        
        for r in results:
            continue
    
    def _download_single_pdb(self, pdbid_and_chains_tuple):
        
        pdbname = pdbid_and_chains_tuple[0]
        chains = pdbid_and_chains_tuple[1]
        
        possible_links = [l.format(pdbname) for l in POSSIBLELINKS]
        
        with self._attempt_download(pdbname):
            response = self._download_data(possible_links)
        
        try:
            downloaded_data = response.read()
        except (AttributeError, UnboundLocalError):  # response is None
            return
        
        pdbdata = PDBDataConstructor(downloaded_data)
        
        pdbdata.build()
        
        if chains[0] is None:
            chains = pdbdata.chain_set
        
        for chain in chains:
            
            pdbdata.add_filter_record_name(self.record_name)
            pdbdata.add_filter_chain(chain)
            destination = Path(self.destination, f'{pdbname}_{chain}.pdb')
            try:
                pdbdata.write(destination)
            except EmptyFilterError:
                log.error(traceback.format_exc())
                log.error(f'Empty Filter for {destination}')
            pdbdata.clear_filters()
        
    def _download_data(self, possible_links):
        for weblink in possible_links:
            try:
                response = urllib.request.urlopen(weblink)
            except urllib.error.HTTPError:
                log.error(S(f'failed from {weblink}'))
                continue
            else:
                log.info(S(f'completed from {weblink}'))
                return response
        else:
            raise DownloadFailedError
    
    @contextlib.contextmanager
    def _attempt_download(self, pdbname):
        try:
            yield
        except DownloadFailedError as e:
            log.error(S(f'{repr(e)}: FAILED {pdbname}'))
            return
