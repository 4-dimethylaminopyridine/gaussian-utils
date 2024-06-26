import argparse
import pathlib
from .conformer_searcher import ConformerSearcher
from .gaussian_input_generator import GaussianInputGenerator
from .gaussian_results_extractor import GaussianResultsExtractor

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

def extract_gs_results_subcommand(args):
    results_directory = pathlib.Path(args.input_dir)
    output_directory = pathlib.Path(args.output_dir)

    if not results_directory.is_dir():
        raise FileNotFoundError(f'{results_directory} does not exist!')

    extractor = GaussianResultsExtractor(results_directory)
    extractor.save_results(output_directory, 'xlsx')

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

    # extract-gs-results subcommand
    subparser_extract_gs_results = subparsers.add_parser(
        'extract-gs-results',
        help='Extract results from Gaussian log files.',
    )
    subparser_extract_gs_results.add_argument(
        'input_dir', metavar='<input_dir>',
        help='path to the directory containing Gaussian log files',
    )
    subparser_extract_gs_results.add_argument(
        'output_dir', metavar='<output_dir>',
        help='path to the output directory',
    )
    subparser_extract_gs_results.set_defaults(func=extract_gs_results_subcommand)


    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
