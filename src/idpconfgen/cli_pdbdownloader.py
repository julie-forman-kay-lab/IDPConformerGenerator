"""
PDBDOWNLOADER
"""
import argparse
from pathlib import Path

ap = argparse.ArgumentParser(
    description=__doc__,
    usage=__doc__,
    )
# https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments 

ap.add_argument(
    'cull',
    help='The CULL PDB file list.',
    type=Path,
    )

def load_args():

    cmd = ap.parse_args()
    return cmd


def main(args):
    print('I am the main of pdbdownloader')
    print(args)

if __name__ == '__main__':
    cmd = load_args()
    main(cmd)
