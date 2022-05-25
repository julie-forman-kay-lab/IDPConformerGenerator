"""
Sample bond geometries from a database.

This database was created from high-resolution PDB structures.
"""
from pathlib import Path


_file = Path(__file__).parent
bgeo_sampling_path = Path(_file, 'bgeo_sampling.tar')
name = "sampling"
