import argparse
import pathlib
from .conformer_searcher import ConformerSearcher
from .gaussian_input_generator import GaussianInputGenerator

def search_conformers_subcommand(args):

    input_file = pathlib.Path(args.input_file)
    dest_dir = pathlib.Path(args.dest_dir)

    # make sure the input file exist
    if not input_file.is_file():
        raise FileNotFoundError(f'{input_file} does not exist!')

    searcher = ConformerSearcher(input_file)
    searcher.generate_results(dest_dir)

def generate_gs_inputs_subcommand(args):

    config_file = pathlib.Path(args.config_file)
    dest_dir = pathlib.Path(args.dest_dir)

    # make sure the config file exist
    if not config_file.is_file():
        raise FileNotFoundError(f'{config_file} does not exist!')

    generator = GaussianInputGenerator(config_file)
    generator.generate_gaussian_input(dest_dir)

def main():

    parser = argparse.ArgumentParser(
        prog='gaussian-utils',
    )
    subparsers = parser.add_subparsers(metavar='subcommands')

    # search-conformers subcommand
    subparser_search_conformers = subparsers.add_parser(
        'search-conformers',
        help='Generate conformers for molecules read from a spreadsheet. The program also produces an html file for users to validate their inputs.'
    )
    subparser_search_conformers.set_defaults(func=search_conformers_subcommand)
    subparser_search_conformers.add_argument(
        'input_file', metavar='<input_file>',
        help='path to the input spreadsheet',
    )
    subparser_search_conformers.add_argument(
        'dest_dir', metavar='<destination_dir>',
        help='path to the output directory',
    )

    # generate-gs-inputs subcommand
    subparser_generate_gs_inputs = subparsers.add_parser(
        'generate-gs-inputs',
        help='Generate Gaussian input files according to the config file.'
    )
    subparser_generate_gs_inputs.add_argument(
        'config_file', metavar='<config_file>',
        help='path to the config file',
    )
    subparser_generate_gs_inputs.add_argument(
        'dest_dir', metavar='<destination_dir>',
        help='path to the output directory',
    )
    subparser_generate_gs_inputs.set_defaults(func=generate_gs_inputs_subcommand)


    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
