import argparse
import pathlib

from .conformer_search import conformer_search
from .molecule_reader import read_from_spreadsheet, spreadsheet_to_html
from .gaussian_input_generator import convert_all


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


def gen_gs_inputs_subcommand(args):
    config_file = pathlib.Path(args.config_file)
    convert_all(config_file)


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

    # gen-gs-inputs subcommand
    subparser_gen_gs_inputs = subparsers.add_parser(
        'gen-gs-inputs',
        help='Generate Gaussian input files according to the config file.'
    )
    subparser_gen_gs_inputs.set_defaults(func=gen_gs_inputs_subcommand)
    subparser_gen_gs_inputs.add_argument(
        'config_file', metavar='<config_file>',
        help='path to the config file'
    )

    args = parser.parse_args()

    args.func(args)


if __name__ == '__main__':
    main()
