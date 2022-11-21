import pathlib
import os

# paths

ROOT_DIR = pathlib.Path(os.path.dirname(
    os.path.realpath(__file__)))
TEMPLATE_DIR = ROOT_DIR.joinpath('templates')

# template name
GAUSSIAN_GJF_TEMPLATE = 'gaussian_input.gjf.j2'

# column names
RDKIT_MOL_COL = 'rdkit_mol'
SMILES_COL = 'SMILES'
NAME_COL = 'Name'
