import pandas
from rdkit.Chem import PandasTools


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
        df, smilesCol='SMILES'
    )
    
    return df