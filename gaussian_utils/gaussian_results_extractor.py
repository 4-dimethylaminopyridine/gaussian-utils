import pathlib
import pandas


class GaussianResultsExtractor:

    def __init__(self, results_directory):

        self._results_directory = pathlib.Path(results_directory)

    def _get_thermo_data(self, gaussian_log_file):

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
            if line.strip().startswith('Sum of electronic and thermal Free Energies'):
                ret['G'] = line.replace(
                    'Sum of electronic and thermal Free Energies=', '').strip()
                ret['G'] = float(ret['G'])

        return ret

    def _hartree_to_kcal_per_mol(self, hartree):
        return hartree * 627.51

    def save_results(self, output_directory, output_format):

        df_dict = {
            'File Name': [],
            'G (hartree)': [],
            'G + 1.89 (kcal/mol)': [],
        }

        for f in self._results_directory.glob('**/*.log'):

            print(f'Extracting {f}')

            file_name = str(f.relative_to(self._results_directory))

            # read thermochemistry data
            thermo_data = self._get_thermo_data(f)
            if 'G' in thermo_data:
                G_hartree = thermo_data['G']
                G_kcal = self._hartree_to_kcal_per_mol(thermo_data['G']) + 1.89
            else:
                G_hartree = None
                G_kcal = None

            df_dict['File Name'].append(file_name)
            df_dict['G (hartree)'].append(G_hartree)
            df_dict['G + 1.89 (kcal/mol)'].append(G_kcal)

        df = pandas.DataFrame(df_dict)
        if output_format == 'xlsx':
            df.to_excel(f'{output_directory}/results.xlsx', index=False)
        else:
            raise ValueError(f'Unsupported output format: {output_format}')
