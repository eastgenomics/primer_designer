'''
Script to design primers and generate a PDF report based off passed
coordiantes using primer3 & smalt.

Original: Kim Brugger (15 Sep 2015), contact: kim@brugger.dk
Made better: Nikita Povarnitsyn (09 Apr 2019)
Made even better: Jethro Rainford (22/03/2021)
'''
import os
import sys
import re
import subprocess
import shlex
import argparse
import reportlab
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

# paths to req. tools, defined in config
from config import SAMTOOLS, TABIX, SMALT, PRIMER3, FONT,\
    REF_37, REF_38, DBSNP_37, DBSNP_38

# Default parameters
FLANK = 500
FUSION = False
TARGET_LEAD = 50
NR_PRIMERS = 4
ALLOWED_MISMATCHES = 0
MAX_MAPPINGS = 5
REFERENCE = None  # depends on the chosen reference genome
DBSNP = None  # depends on the chosen reference genome
VERBOSE = 2
VERSION = '1.2'

TMP_FILES = []

colours = [[255, 0, 0],  # red
           [0, 255, 0],  # green
           [0, 0, 255],  # blue
           [255, 0, 255],  # Pink
           [0, 255, 255],  # cyan
           [255, 255, 0]]  # yellow


class Fusion():
    """Fusion primer designing related functions"""

    def split_input(self, coords):
        """
        Creates nested dict of dicts to store chrom, pos, and the sequence
        sign when designing for fusions

        Args:
            - coords (str): coords of targets
        Returns:
            - coord_dicts (dict): nested dict of coordinate dicts
        """
        coord_dicts = {}
        coordinates = coords.split("_")

        for i, coordinate in enumerate(coordinates):
            instance = coordinate.split(":")

            # replacing the a/b with symbols, cos i wrote the whole script with
            # symbols and was lazy to replace them
            if instance[2] == "a":
                instance[2] = re.sub('a', '<', instance[2])

            elif instance[2] == "b":
                instance[2] = re.sub('b', '>', instance[2])

            coord_dicts[i] = {
                'CHR': instance[0],
                'POS': int(instance[1]),
                'SIDE': instance[2],
                'STRAND': instance[3]}

        return coord_dicts


    def fetch_seqs(self, coord_dict):
        """
        Get the sequences for each coordinates and add the sequences in
        the nested dictionary

        Args:
            - coord_dicts (dict): nested dict of coordinate dicts
        Returns:
            - coord_dicts (dict): nested dict of coordinate dicts
        """
        for index, dict in coord_dict.items():
            if dict['SIDE'] == "<":
                coord_dict[index]['SEQ'] = Sequence.fetch_region(
                    dict['CHR'], dict['POS'], dict['POS'] + FLANK
                )

            elif dict['SIDE'] == ">":
                coord_dict[index]['SEQ'] = Sequence.fetch_region(
                    dict['CHR'], dict['POS'] - FLANK, dict['POS']
                )

            else:
                print(
                    "The fusion sequence was marked incorrectly. \n  \
                    The signs to be used: <>. The sign used:  \
                    dict['SIDE']"
                )
                break

        return coord_dict


    def flip_fusion_seq(self, seqs_dict):
        """
        Flips a sequence in case two same sides relative to breakpoint
        are added

        Args:
            -
        Returns:
            -
        """
        complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            '<': '>', '>': '<', '[': ']', ']': '['
        }

        # DARK_SIDE is requrired to store original SIDE value, which is
        # used in pick_best_primer() as well as printing coordinates in
        # the right order

        for i in range(0, len(seqs_dict)):

            if seqs_dict[i]['SIDE'] == ">" and seqs_dict[i]['STRAND'] == "1":
                seqs_dict[i]['DARK_SIDE'] = ">"

            elif seqs_dict[i]['SIDE'] == "<" and seqs_dict[i]['STRAND'] == "1":
                seqs_dict[i]['DARK_SIDE'] = "<"

            elif seqs_dict[i]['SIDE'] == ">" and seqs_dict[i]['STRAND'] == "-1":
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])
                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])
                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = "<"

            elif seqs_dict[i]['SIDE'] == "<" and seqs_dict[i]['STRAND'] == "-1":
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])
                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])
                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = ">"

        if seqs_dict[0]['DARK_SIDE'] == ">" and seqs_dict[1]['DARK_SIDE'] == "<":
            target_sequence = ''.join([seqs_dict[0]['SEQ'], seqs_dict[1]['SEQ']])
            marked_sequence = ''.join([seqs_dict[0]['MSEQ'], seqs_dict[1]['MSEQ']])
            tagged_string = ''.join([seqs_dict[0]['TSEQ'], seqs_dict[1]['TSEQ']])

        elif seqs_dict[0]['DARK_SIDE'] == "<" and seqs_dict[1]['DARK_SIDE'] == ">":
            target_sequence = ''.join([seqs_dict[1]['SEQ'], seqs_dict[0]['SEQ']])
            marked_sequence = ''.join([seqs_dict[1]['MSEQ'], seqs_dict[0]['MSEQ']])
            tagged_string = ''.join([seqs_dict[1]['TSEQ'], seqs_dict[0]['TSEQ']])

        else:
            print("This fusion is not possible")
            sys.exit()

        return target_sequence, tagged_string, marked_sequence


class Sequence():
    """Sequence related functions"""

    @staticmethod
    def fetch_region(chrom, start, end):
        """
        Gets the nucleotide sequence from the appropriate reference
        file with the given coordinates

        Args:
            - chrom (str): target chromosome
            - start (int): target start coord
            - end (int): target end coord

        Returns:
            - sequence (str): target sequence from reference file

        """
        cmd = f"{SAMTOOLS} faidx {REFERENCE}  {chrom}:{start}-{end} "
        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)

        output = p.communicate()
        sequence = ""

        for line in (output[0].split("\n")):
            if re.match('>', line):
                # skips lines containing '>' for some reason...
                continue

            sequence += line

        return sequence


    def markup_sequence(
        self, flank, sequence, FUSION=False, chrom=None, startpos=None,
        endpos=None
    ):
        """
        Marks the target sequence with "-*-", creates a string of
        characters representing the SNPs, target lead (50 bases up and 
        down the position)

        """
        if FUSION:
            for key, value in sequence.iteritems():
                sequence_list = list(sequence['SEQ'])
                tags = [" "] * len(sequence['SEQ'])
                chrom = sequence['CHR']
                side = sequence['SIDE']

                if side == "<":
                    startpos = int(sequence['POS'])
                    endpos = startpos + FLANK

                    tags[0] = '*'

                    for x in range(0, len(tags)):
                        if x in range(1, 1 + TARGET_LEAD):
                            tags[x] = '-'

                    sequence_list[TARGET_LEAD -1] = sequence_list[TARGET_LEAD - 1] + '] '

                if side == ">":
                    startpos = int(sequence['POS']) - FLANK
                    endpos = int(sequence['POS'])

                    tags[-1] = '*'

                    for x in range(0, len(tags)):
                        if x in range(len(tags) - TARGET_LEAD, len(tags) - 1):
                            tags[x] = '-'

                    sequence_list[len(tags) - TARGET_LEAD] = sequence_list[
                        len(tags) - TARGET_LEAD] + ' ['

                dbSNPs = self.fetch_known_SNPs(DBSNP, chrom, startpos, endpos)
                tags = self.markup_repeats(tags, sequence['SEQ'])

                sequence_list, tagged_string = self.markup_SNPs(
                    dbSNPs, sequence_list, tags, startpos, endpos, FUSION, side
                )

        else:
            # run normal markup sequence
            sequence_list = list(sequence)

            # Creating a list of spaces of size of the fetched sequence
            tags = [" "] * len(sequence)

            # Our target base
            start = FLANK  # FLANK=500 bases
            end = len(tags) - flank

            for x in range(0, len(tags)):
                # Tagging the central 501st base
                if x >= start and x <= end:
                    tags[x] = '*'
                # TARGET_LEAD is 50 bases, 50 bases up nad down the central
                # nucleotide are tagged
                if x in range(start - TARGET_LEAD,
                            start) or x in range(end, end + TARGET_LEAD):
                    tags[x] = '-'

            # Tag the target sequence (50 bases up and down the central
            # nucleotide), regardless of the nature of the variant only one (1)
            # base is tagged as the target.
            sequence_list[flank - TARGET_LEAD] = ' [' + \
                sequence[flank - TARGET_LEAD]

            sequence_list[(- flank + TARGET_LEAD) - 1] = sequence_list[
                (- flank + TARGET_LEAD) - 1] + '] '

            dbSNPs = self.fetch_known_SNPs(
                DBSNP, chrom, startpos - flank, endpos + flank
            )

            tags = self.markup_repeats(tags, sequence)

            sequence_list, tagged_string = self.markup_SNPs(
                dbSNPs, sequence_list, tags, startpos, endpos
            )

        return (sequence_list, tagged_string)


    def markup_repeats(self, tags, sequence):
        """
        marks the repeats in the sequence longer than 6 bases
        """
        regex = r"(.{1,})\1\1\1\1\1+"
        matches = re.finditer(regex, sequence, re.MULTILINE)
        for matchNum, match in enumerate(matches, start=1):

            print("Match {matchNum} was found at {start}-{end}: {match}".format(
                matchNum=matchNum, start=match.start(),
                end=match.end(), match=match.group())
            )

            for x in range(match.start(), match.end()):

                if len(match.group()) >= 10:
                    tags[x] = '^'

                elif len(match.group()) < 10:
                    tags[x] = '+'

        return tags


    def markup_SNPs(self, dbSNPs, sequence_list, tags, startpos,
                    endpos, FUSION=False, side=None):
        """
        Marks the known snps depending on the variables passed
        """

        masked_positions = []

        for dbSNP in dbSNPs:

            if len(dbSNP) < 6:
                continue

            snp_chr, snp_pos, snp_id, snp_ref, snp_alt, common = dbSNP[0:6]

            snp_pos = int(snp_pos)

            if FUSION:
                if side == "<":
                    if common == "0":
                        continue

                    if snp_pos >= startpos + TARGET_LEAD:
                        mask_pos = snp_pos - startpos

                    else:
                        continue

                elif side == ">":
                    if common == "0":
                        continue

                    if snp_pos <= endpos - TARGET_LEAD:
                        mask_pos = snp_pos - startpos

                    else:
                        continue

            else:
                if snp_pos >= startpos - TARGET_LEAD and snp_pos <= endpos + TARGET_LEAD:
                    continue

                if common == '0':
                    continue

                if common == '1':

                    mask_pos = snp_pos - (startpos - FLANK)

            # In the odd case we are looking at common deletion, mask the whole
            # region. Normally this will just be one base
            for i in range(0, len(snp_ref)):

                if len(sequence_list) <= mask_pos + i:
                    break

                # If already masked skip masking it again.
                if re.search('<',sequence_list[mask_pos + i]) or re.search(
                    r'\[',sequence_list[mask_pos + i]):

                    continue

                sequence_list[mask_pos + i] = ' <' + \
                    sequence_list[mask_pos + i] + '> '
                tags[mask_pos + i] = 'X'

            else:
                pass

        sequence_list = "".join(sequence_list)
        sequence_list = re.sub(' ', '', sequence_list)
        tagged_string = "".join(tags)

        return sequence_list, tagged_string


    def fetch_known_SNPs(self, tabix_file, chrom, start, end):
        """
        Getting the known SNPs for the genomic region of interest
        """
        cmd = f"{TABIX} {tabix_file}  {chrom}:{start}-{end}"

        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)

        output = p.communicate()
        var = []

        for line in output[0].split("\n"):
            var.append(line.split("\t"))

        return var




def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chr')
    # only the position or range can be given
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--pos', type=int)
    group.add_argument('-r', '--range', nargs=2)
    group.add_argument(
        '-b', '--blend', type=str,
        help=(
            'The input must be in this format chr1:pos:side:strand_chr2:pos:'
            'side:strand, where side = a or b (after or before breakpoint) '
            'and strand = 1 or -1'
        )
    )

    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank', type=int)
    parser.add_argument(
        '-t', '--text_output', action='store_true',
        help='Saves the report in txt file including the consensus sequence'
    )
    parser.add_argument(
        '--grch37', help='use GRCh37 as reference genome', action="store_true"
    )
    parser.add_argument(
        '--grch38', help='use GRCh38 as reference genome', action="store_true"
    )

    args = parser.parse_args()

    if args.grch37:
        REFERENCE = REF_37
        DBSNP = DBSNP_37

    elif args.grch38:
        REFERENCE = REF_38
        DBSNP = DBSNP_38
    else:
        print("Please select a reference genome to use")
        parser.parse_args('-h')

    return(args, REFERENCE, DBSNP)


def main():

    args, REFERENCE, DBSNP = parse_args()
    fusion = Fusion()
    sequence = Sequence()

    chrom = args.chr

    if args.range:
        startpos, endpos = [int(x) for x in args.range]
        seqs = None

    elif args.blend:
        # blend?
        fusion = args.blend
        FUSION = True
        chrom = None
        startpos = None
        endpos = None

    elif args.pos:
        startpos = args.pos
        endpos = args.pos
        seqs = None

    ### Sequence retrieval and markup ###

    # If fusion is passed, run functions required to get the sequences
    # and mark them
    if FUSION:
        coord_dicts = fusion.split_input(fusion)
        seqs = fusion.fetch_seqs(coord_dicts)

        for ids, dicts in seqs.items():
            marked_seq, tagged_seq = sequence.markup_sequence(FLANK, dicts, FUSION)
            dicts['MSEQ'] = marked_seq
            dicts['TSEQ'] = tagged_seq
            dicts['DARK_SIDE'] = ''

        target_sequence, tagged_string, marked_sequence = fusion.flip_fusion_seq(seqs)
        region_id = fusion
    else:
        # Normal run: when either a position or range is passed
        target_sequence = sequence.fetch_region(chrom, startpos - FLANK, endpos + FLANK)
        marked_sequence, tagged_string = sequence.markup_sequence(
            FLANK, target_sequence, FUSION, chrom, startpos, endpos
        )

        if startpos == endpos:
            region_id = "{}_{}".format(chrom, startpos)
        else:
            region_id = "{}_{}_{}".format(chrom, startpos, endpos)


    # Creating primers and choosing the best ones; and converting to strings
    primer3_results = run_primer3(region_id, marked_sequence)

    passed_primers = check_primers(
        region_id, target_sequence, primer3_results, chrom,
        startpos, endpos, FUSION, seqs
    )

    passed_primer_seqs = extract_passed_primer_seqs(
        primer3_results, passed_primers
    )

    if not FUSION:
        # Converting the primers into strings which can be out into the
        # PDF report
        mapped_primer_strings, mapped_primer_colours = make_primer_mapped_strings(
            target_sequence, passed_primer_seqs
        )

        lines = pretty_print_mappings(
            target_sequence, tagged_string, mapped_primer_strings,
            startpos - FLANK
        )

    # generate output PDF report
    filename = re.sub("[<>:]", "_", region_id)

    if args.output:
        filename = args.output 

    if args.text_output:
        pretty_primer_data(
            "{}.txt".format(
                filename, primer3_results, passed_primers, chrom, startpos,
                endpos, fwd_primer, rev_primer, target_sequence, FUSION, seqs
            )
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
        top_offset = pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width,
            FUSION, chrom, startpos, endpos, seqs
        )

        c.setFont('mono', 8)
        pretty_pdf_fusion_mappings(
            top_offset, c, seqs, passed_primer_seqs, FUSION)
        pretty_pdf_method(top_offset, args, c, seqs)

    else:
        # normal primer report formatting
        c.setFont('mono', 6)
        top_offset = pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width, FUSION,
            chrom, startpos, endpos
        )
        c.setFont('mono', 8)
        pretty_pdf_mappings(
            top_offset, target_sequence, tagged_string, mapped_primer_strings,
            mapped_primer_colours, startpos - FLANK, c
        )
        pretty_pdf_method(top_offset, args, c)

    c.showPage()
    c.save()

    for filename in TMP_FILES:
        print("deleting tmp file: {}".format(filename))
        os.remove(filename)


if __name__ == '__main__':
    main()
