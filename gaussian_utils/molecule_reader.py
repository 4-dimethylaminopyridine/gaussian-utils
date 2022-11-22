import pathlib
import pandas
from rdkit.Chem import PandasTools

from gaussian_utils.config import *


def read_from_spreadsheet(path, spreadsheet_format):
    """Read a list of molecules from a spreadsheet.

    Parameters
    ----------
    path : pathlib.Path
        Path to the spreadsheet file.
    spreadsheet_format: str
        csv or excel.

    Returns
    -------
    A pandas dataframe containing a rdkit col.

    """
    if spreadsheet_format == 'excel':
        df = pandas.read_excel(path)
    elif spreadsheet_format == 'csv':
        df = pandas.read_csv(path)
    else:
        raise ValueError('Wrong spreadsheet_format "%s"!' % spreadsheet_format)

    PandasTools.AddMoleculeColumnToFrame(
        df, smilesCol=SMILES_COL, molCol=RDKIT_MOL_COL
    )

    return df


def spreadsheet_to_html(path, spreadsheet_format):

    output = path.parent.joinpath(path.stem + '.html')
    df = read_from_spreadsheet(path, spreadsheet_format)
    df.to_html(str(output))
