"""Test the side-chain packing components."""
import sys
import types

import numpy as np

from idpconfgen.components.sidechain_packing import init_mcsce_sidechains


def _install_fake_mcsce(monkeypatch):
    """
    Install a stand-in for the MCSCE package.

    Like MCSCE, `initialize_func_calc` keeps the state of one sequence in
    module globals, which `create_side_chain` reads. Its structures have
    one side-chain atom (CB) per residue.
    """
    builder = types.ModuleType('mcsce.core.side_chain_builder')
    builder.sidechain_placeholders = []
    builder.energy_calculators = []

    def initialize_func_calc(efunc_creator, structure=None):
        builder.sidechain_placeholders = [structure.fasta]
        builder.energy_calculators = [efunc_creator]

    def create_side_chain(structure, **params):
        names = ['N', 'CA', 'C', 'O', 'CB'] * len(builder.sidechain_placeholders[0])  # noqa: E501
        return types.SimpleNamespace(
            data_array=np.array([[0, 0, name] for name in names]),
            coords=np.arange(3 * len(names), dtype=float).reshape(-1, 3),
            )

    builder.initialize_func_calc = initialize_func_calc
    builder.create_side_chain = create_side_chain

    class Structure:
        def __init__(self, fasta):
            self.fasta = fasta

        def build(self):
            pass

    definitions = types.ModuleType('mcsce.core.build_definitions')
    definitions.forcefields = {'Amberff14SB': lambda **kwargs: None}
    libenergy = types.ModuleType('mcsce.libs.libenergy')
    libenergy.prepare_energy_function = lambda **kwargs: None
    libstructure = types.ModuleType('mcsce.libs.libstructure')
    libstructure.Structure = Structure
    core = types.ModuleType('mcsce.core')
    core.build_definitions = definitions
    core.side_chain_builder = builder
    libs = types.ModuleType('mcsce.libs')
    libs.libenergy = libenergy
    libs.libstructure = libstructure
    package = types.ModuleType('mcsce')
    package.core = core
    package.libs = libs
    modules = {
        'mcsce': package,
        'mcsce.core': core,
        'mcsce.core.build_definitions': definitions,
        'mcsce.core.side_chain_builder': builder,
        'mcsce.libs': libs,
        'mcsce.libs.libenergy': libenergy,
        'mcsce.libs.libstructure': libstructure,
        }
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)


def test_mcsce_sequences_in_one_process(monkeypatch):
    """Test each sequence packs with its own MCSCE state and atom mask."""
    _install_fake_mcsce(monkeypatch)
    masks = types.SimpleNamespace(non_sidechains=slice(None), all_sidechains='mask')  # noqa: E501
    short = init_mcsce_sidechains('AAA', masks, masks)
    long = init_mcsce_sidechains('AAAAA', masks, masks)
    coords = np.zeros((10, 3))
    # the side chains of the first sequence after the second initialized
    # MCSCE, then of the second one
    for calc, residues in ((short, 3), (long, 5), (short, 3), (long, 5)):
        mask, sidechains = calc(coords)
        assert mask == 'mask'
        assert sidechains.shape == (residues, 3)
