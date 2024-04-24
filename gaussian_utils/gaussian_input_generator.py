import pathlib
import yaml
import rdkit.Chem
import jinja2


class GaussianInputGenerator:

    _GAUSSIAN_GJF_TEMPLATE = '''{{ gaussian_route_section }}

{{ gaussian_title_section }}

{{ molecule_charge }} {{ molecule_spin }}
{{ molecule_xyz }}
{% if optional_additional_sections is not none %}{{ optional_additional_sections }}{% endif %}
! This is a comment line to indicate a necessary blank line
'''

    def __init__(self, config_file_path):

        self._config_file_path = pathlib.Path(config_file_path)

    def _sdf_to_gjf(self, sdf_path, output_dir, additional_args):

        """Generate a list of gjf files from a sdf file. Each entry in the sdf file
        must contain a "ID" property.

        Parameters
        ----------
        sdf_path : pathlib.Path
            Path to the sdf file.
        output_dir : pathlib.Path
            Path to the output directory.
        additional_args: dict
            additional arguments required for template rendering.

        Returns
        -------
        None.

        """

        # create directory
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(parents=True)

        reader = rdkit.Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        for mol in reader:
            ID = mol.GetProp('ID')
            molecule_xyz = rdkit.Chem.MolToXYZBlock(mol)
            molecule_xyz = ''.join(molecule_xyz.splitlines(keepends=True)[2:]) # need to remove the first two lines
            molecule_xyz = molecule_xyz.strip()

            args = additional_args.copy()
            args['molecule_xyz'] = molecule_xyz

            self._render_template(
                self._GAUSSIAN_GJF_TEMPLATE,
                output_dir / f'ID_{ID}.gjf',
                args,
            )

    def mol2_to_gjf(self, mol2_path, output_file, additional_args):

        mol = rdkit.Chem.MolFromMol2File(str(mol2_file), removeHs=False)

        molecule_xyz = rdkit.ChemMolToXYZBlock(mol)

        # need to remove first 2 lines
        molecule_xyz = ''.join(molecule_xyz.splitlines(keepends=True)[2:])
        molecule_xyz = molecule_xyz.strip()

        args = additional_args.copy()
        args['molecule_xyz'] = molecule_xyz

        self._render_template(
            self._GAUSSIAN_GJF_TEMPLATE,
            output_file,
            args
        )

    def _render_template(self, template_str, dest_file, args):
        jinja2_template = jinja2.Environment(
            loader=jinja2.BaseLoader()
        ).from_string(template_str)

        with open(dest_file, 'w') as f:
            f.write(jinja2_template.render(**args))

    def generate_gaussian_input(self, output_dir):

        output_dir = pathlib.Path(output_dir)

        with self._config_file_path.open() as f:
            yaml_data = yaml.safe_load(f)

        for item in yaml_data['files']:

            filename = item['filename']
            output_name = item['output_name']
            gaussian_title_section = item['gaussian_title_section']
            gaussian_route_section = item['gaussian_route_section']
            molecule_charge = item['molecule_charge']
            molecule_spin = item['molecule_spin']

            # can leave this blank
            if 'optional_additional_sections' in item:
                optional_additional_sections = [i.strip() for i in item['optional_additional_sections']]
                optional_additional_sections = '\n' + '\n\n'.join(optional_additional_sections)
            else:
                optional_additional_sections = None

            template_args =  {
                'gaussian_title_section': gaussian_title_section,
                'gaussian_route_section': gaussian_route_section,
                'molecule_charge': molecule_charge,
                'molecule_spin': molecule_spin,
                'optional_additional_sections': optional_additional_sections
            }

            for f in self._config_file_path.parent.glob(filename):
                print(f'Converting {f}')
                if f.suffix == '.sdf':
                    self._sdf_to_gjf(f, f'{output_dir}/{output_name}', template_args)
                elif f.suffix == '.mol2':
                    self._mol2_to_gjf(f, f'{output_dir}/{output_name}.gjf', template_args)
                else:
                    raise ValueError(f'Unsupported file type: {f.suffix}')
