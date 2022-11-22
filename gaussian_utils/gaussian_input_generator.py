import yaml
import jinja2
import rdkit.Chem

from .config import *


def render_template(template, dest_file, args):
    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(TEMPLATE_DIR))
    )
    template = jinja_env.get_template(template)

    f = open(dest_file, 'w')
    f.write(template.render(**args))
    f.close()


def sdf_to_gjf(sdf_path, output_dir, additional_args):
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
    output_dir.mkdir(parents=True, exist_ok=True)

    reader = rdkit.Chem.SDMolSupplier(str(sdf_path))
    for mol in reader:
        ID = mol.GetProp('ID')
        molecule_xyz = rdkit.Chem.MolToXYZBlock(mol)
        # need to remove first 2 lines
        molecule_xyz = ''.join(molecule_xyz.splitlines(keepends=True)[2:])

        args = additional_args.copy()
        args['molecule_xyz'] = molecule_xyz

        render_template(
            GAUSSIAN_GJF_TEMPLATE,
            output_dir.joinpath('ID_%s.gjf' % ID),
            args
        )


def mol2_to_gjf(mol2_path, output_file, additional_args):
    mol = rdkit.Chem.MolFromMol2File(str(mol2_path), removeHs=False)

    molecule_xyz = rdkit.Chem.MolToXYZBlock(mol)
    # need to remove first 2 lines
    molecule_xyz = ''.join(molecule_xyz.splitlines(keepends=True)[2:])

    args = additional_args.copy()
    args['molecule_xyz'] = molecule_xyz

    render_template(
        GAUSSIAN_GJF_TEMPLATE,
        output_file,
        args
    )


def convert_all(config_file_path):

    with open(str(config_file_path), 'r') as f:
        yaml_data = yaml.safe_load(f)

    for item in yaml_data['files']:

        filename = item['filename']
        gaussian_title_section = item['gaussian_title_section']
        gaussian_route_section = item['gaussian_route_section']
        molecule_charge = item['molecule_charge']
        molecule_spin = item['molecule_spin']

        template_args = {
            'gaussian_title_section': gaussian_title_section,
            'gaussian_route_section': gaussian_route_section,
            'molecule_charge': molecule_charge,
            'molecule_spin': molecule_spin
        }

        for f in config_file_path.parent.glob(filename):
            print('Converting %s' % f)
            if f.suffix == '.sdf':
                output_dir = f.parent.joinpath(f.stem + '_gjf')
                sdf_to_gjf(f, output_dir, template_args)
            elif f.suffix == '.mol2':
                output_file = f.parent.joinpath(f.stem + '.gjf')
                mol2_to_gjf(f, output_file, template_args)
            else:
                raise ValueError(
                    'Unsupported file type ("%s") matched!' % f.name)
