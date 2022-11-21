import yaml
import jinja2
import rdkit.Chem

from gaussian_utils.config import *


def render_template(template, dest_file, args):
    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(TEMPLATE_DIR))
    )
    template = jinja_env.get_template(template)

    f = open(dest_file, 'w')
    f.write(template.render(**args))
    f.close()


def convert_all(config_file_path):

    def sdf_to_gjf(sdf_path, output_dir):
        # create directory
        output_dir.mkdir(parents=True, exist_ok=True)

        reader = rdkit.Chem.SDMolSupplier(str(sdf_path))
        for mol in reader:
            conformer_ID = mol.GetProp('ID')
            molecule_xyz = rdkit.Chem.MolToXYZBlock(mol)
            # need to remove first 2 lines
            molecule_xyz = ''.join(molecule_xyz.splitlines(keepends=True)[2:])

            args = {
                'gaussian_title_section': gaussian_title_section,
                'gaussian_route_section': gaussian_route_section,
                'molecule_charge': molecule_charge,
                'molecule_spin': molecule_spin,
                'molecule_xyz': molecule_xyz
            }

            render_template(
                GAUSSIAN_GJF_TEMPLATE,
                output_dir.joinpath('conformer_%s.gjf' % conformer_ID),
                args
            )

    with open(str(config_file_path), 'r') as f:
        yaml_data = yaml.safe_load(f)

    for item in yaml_data['files']:

        filename = item['filename']
        gaussian_title_section = item['gaussian_title_section']
        gaussian_route_section = item['gaussian_route_section']
        molecule_charge = item['molecule_charge']
        molecule_spin = item['molecule_spin']

        for f in config_file_path.parent.glob(filename):
            if f.suffix == '.sdf':
                output_dir = f.parent.joinpath(f.stem)
                sdf_to_gjf(f, output_dir)
