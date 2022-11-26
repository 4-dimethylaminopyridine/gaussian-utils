import cclib
import rdkit.Chem
import pandas

from gaussian_utils.config import *


def get_thermo_data(gaussian_log_file):
    """Extract thermochemistry data from a gaussian log file.

    Parameters
    ----------
    gaussian_log_file : pathlib.Path
        Path to the gaussian log file.

    Returns
    -------
    dict
        A dict containing thermochemistry data.

    """

    ret = {}

    for line in open(str(gaussian_log_file), 'r').readlines():
        if line.strip().startswith('Sum of electronic and thermal Enthalpies'):
            ret['H'] = line.replace(
                'Sum of electronic and thermal Enthalpies=', '').strip()
        if line.strip().startswith('Sum of electronic and thermal Free Energies'):
            ret['G'] = line.replace(
                'Sum of electronic and thermal Free Energies=', '').strip()

    return ret


# TODO
def get_mol_3d(gaussian_log_file):

    tmp_file_path = '/tmp/gaussian_utils.tmp'
    tmp_file = open(tmp_file_path, 'w')

    # convert gaussian log file into xyz
    parser = cclib.io.ccopen(gaussian_log_file)
    data = parser.parse()
    cclib.io.ccwrite(data, 'xyz', tmp_file)
    tmp_file.close()

    # parse tmp xyz file
    mol = rdkit.Chem.MolFromXYZFile(tmp_file_path)

    return {RDKIT_MOL_COL: mol}


def save_results(directory, output_format):

    # hartree to kcal/mol conversion
    # TODO: whether to use 1.89 correction still requires some justification
    def hartree_to_kcal_per_mol(hartree):
        return hartree * 627.51 + 1.89

    dest_html = directory.joinpath('results.html')
    dest_excel = directory.joinpath('results.xlsx')

    df = pandas.DataFrame(
        columns=['File name', 'H', 'G']
    )

    for f in directory.glob(r'**/*.log'):

        print('Extracting %s' % str(f))

        results = {
            'File name': str(f.relative_to(directory))
        }

        # read thermochemistry data
        thermo_data = get_thermo_data(f)
        for key, value in thermo_data.items():
            thermo_data[key] = hartree_to_kcal_per_mol(float(value))
        results.update(thermo_data)

        # write to df
        results_df = pandas.DataFrame(results, index=[0])
        df = pandas.concat([df, results_df])

    if output_format == 'html':
        df.to_html(str(dest_html), index=False)
    elif output_format == 'excel':
        df.to_excel(str(dest_excel), index=False)
