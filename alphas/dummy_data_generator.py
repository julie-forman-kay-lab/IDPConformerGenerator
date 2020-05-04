import random
import string
import sys

aa3to1 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

aminoacids = list(aa3to1.values())

SSletters = ['L', 'H', 'S']


class NameGenerator:
    def __init__(self):
        self.names = []

    def generate_fout_name(self):
        name = self._generate()
        while name in self.names:
            name = self._generate()
        
        self.names.append(name)
        return name
          
    def _generate(self):
        letters = [random.choice(string.ascii_uppercase) for i in range(4)]
        return '{}{}{}{}_{}.data'.format(random.randint(0, 9), *letters)


def generate_fasta():
    protein_length = random.randint(30, 400)
    return [random.choice(aminoacids) for i in range(protein_length)]


def generate_secondary_structure(length):
    
    ss = []
    while len(ss) < length:
        loop_len = random.randint(2, 25)
        ss.extend([random.choice(SSletters)] * loop_len)
    else:
        return ss[:length]


def generate_angles(length):
    
    angles = []
    for residue in range(length):
    
        angles.append(','.join([
            '{:.3f}'.format(random.uniform(-180, 180)),
            '{:.3f}'.format(random.uniform(-180, 180)),
            '{:.3f}'.format(random.uniform(-180, 180)),
            '{:.3f}'.format(random.uniform(-180, 180)),
            ])
            )

    return angles
        
        
def generate_CA_coords(length):
    
    coords = []
    for residue in range(length):
        coords.append(','.join([
            '{:.3f}'.format(random.uniform(0, 100)),
            '{:.3f}'.format(random.uniform(0, 100)),
            '{:.3f}'.format(random.uniform(0, 100)),
            ]))

    return coords


def generate_dummy_data():
    
    name_generator = NameGenerator()

    fout_name = name_generator.generate_fout_name()

    fasta = generate_fasta()
    ss = generate_secondary_structure(len(fasta))
    angles = generate_angles(len(fasta))
    ca_xyz = generate_CA_coords(len(fasta))
    
    data = (
        fasta,
        ss,
        angles,
        ca_xyz,
        )
    
    assert all(isinstance(l, list) for l in data)

    for i in data:
        assert len(i) == len(fasta), i


    data_residues = [','.join(residue) for residue in zip(*data)]
    data_file = '\n'.join(data_residues)

    with open(fout_name, 'w') as fout:
        fout.write(data_file)


def main(num_pdbs):
    for i in range(num_pdbs):
        generate_dummy_data()


if __name__ == "__main__":
    num_pdbs = int(sys.argv[1])
    main(num_pdbs)
