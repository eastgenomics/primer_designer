'''
Script to design primers and generate a PDF report based off passed
coordiantes using primer3 & smalt.

Original: Kim Brugger (15 Sep 2015), contact: kim@brugger.dk
Made better: Nikita Povarnitsyn (09 Apr 2019)
Made even better: Jethro Rainford (23 Mar 2021)
'''
import argparse
from configparser import ConfigParser
import os
from pathlib import Path
import re
from shutil import which

from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

from fusion import Fusion
from primer import Primer3
from report import Report
from sequence import Sequence
from utils import fetch_region

from constants import FLANK

# font file in static dir
FONT = TTFont(
    'mono', (
        f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}'
        '/static/LiberationMono-Regular.ttf'
    )
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--config', required=False,
        help=(
            'Config file with paths to required reference files, use env '
            'variables if not used.'
        )
    )
    chr_arg = parser.add_argument('-c', '--chr')

    # only the position or range can be given
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-p', '--pos', type=int)
    group.add_argument('-r', '--range', nargs=2)

    fusion_arg = parser.add_argument(
        '--fusion', action="store_true",
        help=(
            "Design primers around breakpoint. --b1 and --b2 must be specified"
            " and in this format chr:pos:side:strand, where side = a or b "
            "(after or before breakpoint) and strand = 1 or -1"
        )
    )
    parser.add_argument(
        '--b1', help="first region to design primers for"
    )
    parser.add_argument(
        '--b2', help="second region to design primers for"
    )

    parser.add_argument('-o', '--output', help='output filename prefix')
    parser.add_argument('-d', '--directory', help='directory to generate pdfs')
    parser.add_argument(
        '-f', '--flank', type=int,
        help='flank region to design within (default: 500)'
    )
    parser.add_argument(
        '-t', '--text_output', action='store_true',
        help='Saves the report in txt file including the consensus sequence'
    )
    ref_arg = parser.add_argument(
        '--grch37', help='use GRCh37 as reference genome', action="store_true"
    )
    parser.add_argument(
        '--grch38', help='use GRCh38 as reference genome', action="store_true"
    )

    args = parser.parse_args()

    if not args.grch37 and not args.grch38:
        parser.print_help()
        print('')

        raise argparse.ArgumentError(
            ref_arg, 'Please select a reference genome to use'
        )

    if args.chr and not (args.pos or args.range):
        parser.print_help()
        print('')

        raise argparse.ArgumentError(
            chr_arg, 'chromsome given with no position or range'
        )

    if args.fusion and not (args.b1 and args.b2):
        parser.print_help()
        print('')

        raise argparse.ArgumentError(
            fusion_arg, 'both --b1 AND --b2 must be given for fusion designs.'
        )

    if args.fusion:
        # check b1 and b2 in correct format
        format_regex = '[A-Za-z0-9]*:[0-9]*:[ab]:[-]?[0-9]'
        b1_match = re.match(format_regex, args.b1)
        b2_match = re.match(format_regex, args.b2)

        if b1_match and b2_match:
            b1_match = b1_match.group(0)
            b2_match = b2_match.group(0)

        assert args.b1 == b1_match and args.b2 == b2_match, (
            'ERROR: b1 and / or are in the wrong format. Expected '
            f'chr:pos:side:strand. b1 = {args.b1}. b2 = {args.b2}'
        )

    return args


def load_config(args):
    """
    Load in config parameters, takes either config from args or env variables

    Args:
        - args: cmd line arguments
    Returns:
        - REFERENCE (str): path to appropriate passed reference build file
        - DBSNP (str): path to appropriate passed dbsnp file
    """
    if args.config:
        # passed config file
        config = ConfigParser()
        config.read(args.config)

        if args.grch37:
            REFERENCE = config['REFERENCE']['REF_37']
            SNP = config['REFERENCE']['SNP_37']
            raw_txt = config['REFERENCE']['PRIMER37_TEXT']
        else:
            REFERENCE = config['REFERENCE']['REF_38']
            SNP = config['REFERENCE']['SNP_38']
            raw_txt = config['REFERENCE']['PRIMER38_TEXT']
    else:
        # no config file, Docker
        if args.grch37:
            try:
                REFERENCE = os.environ['REF_37']
                SNP = os.environ['SNP_37']
                raw_txt = os.environ['PRIMER37_TEXT']
            except Exception as e:
                raise ValueError(f'Missing env {e}')
        else:
            REFERENCE = os.environ['REF_38']
            SNP = os.environ['SNP_38']
            raw_txt = os.environ['PRIMER38_TEXT']
    try:
        VERSION = os.environ.get('PRIMER_VERSION', False)
        PRIMER_TEXT = get_chunks(raw_txt, 100)
    except Exception as e:
        raise ValueError(f'Missing env {e}')

    return REFERENCE, SNP, VERSION, PRIMER_TEXT


def get_chunks(sentence: str, maxlength: int):
    """
    Credit: https://stackoverflow.com/questions/57023348/python-splitting-
        a-long-text-into-chunks-of-strings-given-character-limit
    Function to chunk long sentence into smaller chunk

    Input:
        sentence: input sentence
        maxlength: chunk length

    Return:
        list object generator

    """
    start = 0
    end = 0
    while start + maxlength < len(sentence) and end != -1:
        end = sentence.rfind(" ", start, start + maxlength + 1)
        yield sentence[start:end]
        start = end + 1
    yield sentence[start:]


def main():
    # parse args, load in file paths from config
    args = parse_args()
    REFERENCE, SNP, VERSION, PRIMER_TEXT = load_config(args)
    FUSION = False

    # check required tools installed and on PATH
    for tool in ['samtools', 'tabix', 'primer3_core', 'smalt']:
        assert which(tool), f'{tool} is not on path, is it installed?'

    chrom = args.chr

    if args.range:
        startpos, endpos = [int(x) for x in args.range]
        seqs = None

    elif args.fusion:
        fusion_coords = f'{args.b1}_{args.b2}'
        FUSION = True
        chrom = None
        startpos = None
        endpos = None

    elif args.pos:
        startpos = args.pos
        endpos = args.pos
        seqs = None

    # initialise classes
    fusion = Fusion(REFERENCE)
    sequence = Sequence(REFERENCE, SNP)
    primer3 = Primer3(REFERENCE)
    report = Report(FUSION, REFERENCE, PRIMER_TEXT, VERSION, SNP)

    tmp_files = []

    # Sequence retrieval and markup
    if FUSION:
        # If fusion is passed, run functions required to get the sequences
        # and mark them
        coord_dicts = fusion.split_input(fusion_coords)
        seqs = fusion.fetch_seqs(coord_dicts)

        for _, dicts in seqs.items():
            marked_seq, tagged_seq = sequence.markup_sequence(
                FLANK, dicts, FUSION
            )
            dicts['MSEQ'] = marked_seq
            dicts['TSEQ'] = tagged_seq
            dicts['DARK_SIDE'] = ''

        target_sequence, tagged_string, marked_sequence = \
            fusion.flip_fusion_seq(seqs)

        region_id = fusion_coords
    else:
        # Normal run: when either a position or range is passed
        target_sequence = fetch_region(
            chrom,
            startpos - FLANK,
            endpos + FLANK,
            REFERENCE
        )
        marked_sequence, tagged_string = sequence.markup_sequence(
            FLANK, target_sequence, FUSION, chrom, startpos, endpos
        )

        if startpos == endpos:
            region_id = "{}_{}".format(chrom, startpos)
        else:
            region_id = "{}_{}_{}".format(chrom, startpos, endpos)

    # Creating primers, choosing the best ones and converting to strings
    primer3_results = primer3.run_primer3(
        region_id, marked_sequence, tmp_files)

    passed_primers = primer3.check_primers(
        region_id, target_sequence, primer3_results, chrom,
        startpos, endpos, tmp_files, FUSION, seqs
    )

    passed_primer_seqs = primer3.extract_passed_primer_seqs(
        primer3_results, passed_primers
    )

    if not FUSION:
        # Converting the primers into strings which can be out into the
        # PDF report
        mapped_primer_strings, mapped_primer_colours = \
            report.make_primer_mapped_str(
                target_sequence, passed_primer_seqs
            )

    # generate output PDF report, write to output dir in folder
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if args.directory:
        output_dir = f'{parent_dir}/output/{args.directory}/'
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    else:
        output_dir = f'{parent_dir}/output/'

    if args.output:
        filename = f'{output_dir}{args.output}'
    else:
        filename = f'{output_dir}{re.sub("[<>:]", "_", region_id)}'

    if args.text_output:
        filename = f'{output_dir}{filename}.txt'
        report.pretty_primer_data(
            filename, primer3_results, passed_primers, chrom, startpos,
            endpos, target_sequence, FUSION, seqs
        )
    else:
        filename = filename + ".pdf"

    # set PDF formatting
    c = canvas.Canvas(filename, pagesize=A4)
    width, height = A4
    pdfmetrics.registerFont(FONT)

    if FUSION:
        # fusion report, pass required formatting
        c.setFont('mono', 7)
        top_offset = report.pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width, args,
            FUSION, chrom, startpos, endpos, seqs
        )

        c.setFont('mono', 8)
        report.pretty_pdf_fusion_mappings(
            top_offset, c, seqs, passed_primer_seqs, FUSION
        )

        report.pretty_pdf_method(top_offset, args, c, seqs)

    else:
        # normal primer report formatting
        c.setFont('mono', 6)
        top_offset = report.pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width, args,
            FUSION, chrom, startpos, endpos
        )
        c.setFont('mono', 8)
        report.pretty_pdf_mappings(
            top_offset, target_sequence, tagged_string, mapped_primer_strings,
            mapped_primer_colours, startpos - FLANK, c
        )
        report.pretty_pdf_method(top_offset, args, c)

    c.showPage()
    c.save()

    for filename in tmp_files:
        # clean up temporary files
        print("deleting tmp file: {}".format(filename))
        os.remove(filename)


if __name__ == '__main__':
    main()
