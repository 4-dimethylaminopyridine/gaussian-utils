import argparse
import pathlib

from .conformer_search import conformer_search
from .molecule_reader import read_from_spreadsheet, spreadsheet_to_html


def gen_conformers_subcommand(args):

    input_file = pathlib.Path(args.input_file)
    dest_dir = pathlib.Path(args.dest_dir)

    if input_file.suffix == '.xlsx':
        spreadsheet_format = 'excel'
    elif input_file.suffix == '.csv':
        spreadsheet_format = 'csv'
    else:
        # let the read_from_spreadsheet method raise the exception
        spreadsheet_format = input_file.suffix

    df = read_from_spreadsheet(input_file, spreadsheet_format)
    spreadsheet_to_html(input_file, spreadsheet_format)
    conformer_search(df, dest_dir)

def main():

    parser = argparse.ArgumentParser(
        prog='gaussian-utils'
    )

    subparsers = parser.add_subparsers(metavar='subcommands')

    # gen-conformers subcommand
    subparser_gen_conformers = subparsers.add_parser(
        'gen-conformers',
        help='Generate conformers for molecules read from a spreadsheet. The program also produces an html file for users to validate their spreadsheets.'
    )
    subparser_gen_conformers.set_defaults(func=gen_conformers_subcommand)
    subparser_gen_conformers.add_argument(
        'input_file', metavar='<input_file>',
        help='path to the spreadsheet file'
    )
    subparser_gen_conformers.add_argument(
        'dest_dir', metavar='<destination_dir>',
        help='path to the output directory'
    )

    args = parser.parse_args()

    args.func(args)


if __name__ == '__main__':
    main()
