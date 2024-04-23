import pathlib
import pandas
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Draw


class ConformerSearcher:

    _INPUT_NAME_COL = 'Name'
    _INPUT_TYPE_COL = 'Input Type'
    _INPUT_VALUE_COL = 'Input Value'
    _RDKIT_MOL_COL = 'Rdkit Mol'

    def __init__(self, input_path, search_params={}):

        self._input_path = pathlib.Path(input_path)

        self._search_params = {
            'ETKDG_max_iterations': 10000,
            'ETKDG_prune_rms_thresh': 0.5,
            'num_confs': 100
        }
        self._search_params.update(search_params)

        if self._input_path.suffix == '.csv':
            self._df = pandas.read_csv(self._input_path)
        elif self._input_path.suffix == '.xlsx':
            self._df = pandas.read_excel(self._input_path)
        else:
            raise ValueError(f'Unsupported input format: {self._input_path.suffix}')


        # TODO: check for duplicate name

    def _generate_rdkit_mol(self):

        self._df[self._RDKIT_MOL_COL] = None

        for i in self._df.index:

            input_type = self._df[self._INPUT_TYPE_COL][i]
            input_value = self._df[self._INPUT_VALUE_COL][i]

            if input_type == 'SMILES':
                rdkit_mol = rdkit.Chem.MolFromSmiles(input_value)
            else:
                raise ValueError(f'Unsupported input type: "{input_type}"')

            self._df.loc[i, self._RDKIT_MOL_COL] = rdkit_mol

    def _search_conformers(self, m):

        '''Conformer search with ETKDG and then MMFF94.'''

        params = rdkit.Chem.AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.useSmallRingTorsions = True
        params.maxIterations = self._search_params['ETKDG_max_iterations']
        params.pruneRmsThresh = self._search_params['ETKDG_prune_rms_thresh']

        mol = rdkit.Chem.AddHs(m)

        rdkit.Chem.AllChem.EmbedMultipleConfs(
            mol,
            numConfs=self._search_params['num_confs'],
            params=params
        )

        rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs(mol)

        return mol

    def _write_conformers_to_sdf(self, mol, dest):

        writer = rdkit.Chem.SDWriter(str(dest))

        id = 1
        for conf in mol.GetConformers():
            mol.SetProp('ID', str(id))
            writer.write(mol, confId=conf.GetId())
            id += 1

    def _rdkit_mol_to_svg_str(self, rdkit_mol):

        mol = rdkit.Chem.Mol(rdkit_mol)
        mol = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol)
        drawer = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(-1, -1)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_str = drawer.GetDrawingText()
        return svg_str.replace('\n', '')

    def generate_results(self, dest_dir):

        # create dest dir
        dest_dir = pathlib.Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)

        self._generate_rdkit_mol()

        # generate conformers as sdf files
        for i in self._df.index:

            name = self._df[self._INPUT_NAME_COL][i]
            m = self._df[self._RDKIT_MOL_COL][i]
            m_with_conformers = self._search_conformers(m)
            self._write_conformers_to_sdf(
                m_with_conformers,
                f'{dest_dir}/{name}.sdf',
            )

        # generate a reference html file
        self._df[self._RDKIT_MOL_COL] = self._df[self._RDKIT_MOL_COL].apply(self._rdkit_mol_to_svg_str)
        with (dest_dir / 'input.html').open('w') as f:
            f.write(self._df.to_html(escape=False))
