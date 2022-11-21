import rdkit.Chem
import rdkit.Chem.AllChem

from gaussian_utils.config import *


def generate_conformers(mol, numConfs=100):
    """Generate conformers with ETKDG and then MMFF94.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule to do the conformer search.
    numConfs : int, optional
        Maximum number of produced conformers. The default is 100.

    Returns
    -------
    results : rdkit.Chem.rdchem.Mol
        A rdkit Mol object with conformers.

    """

    params = rdkit.Chem.AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.useSmallRingTorsions = True
    params.maxIterations = 10000
    params.pruneRmsThresh = 0.5

    m_with_H = rdkit.Chem.AddHs(mol)
    rdkit.Chem.AllChem.EmbedMultipleConfs(
        m_with_H,
        numConfs=numConfs,
        params=params
    )

    rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs(m_with_H)
    return m_with_H


def write_conformers_to_sdf(mol, dest):

    writer = rdkit.Chem.SDWriter(str(dest))

    id = 1
    for conf in mol.GetConformers():
        writer.write(mol, confId=conf.GetId())


def conformer_search(df, dest_dir):

    # create dest dir
    dest_dir.mkdir(parents=True, exist_ok=True)

    for i in df.index:
        row = df.loc[i]
        name = row[NAME_COL]
        smiles = row[SMILES_COL]

        # generate a new Mol instead of modifying the RDKIT_MOL_COL
        mol = rdkit.Chem.MolFromSmiles(smiles)

        m_with_conformers = generate_conformers(mol)
        write_conformers_to_sdf(
            m_with_conformers, dest_dir.joinpath(name + '.sdf'))
