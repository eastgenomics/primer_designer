'''
Script to design primers and generate a PDF report based off passed
coordiantes using primer3 & smalt.

Original: Kim Brugger (15 Sep 2015), contact: kim@brugger.dk
Made better: Nikita Povarnitsyn (09 Apr 2019)
Made even better: Jethro Rainford (23 Mar 2021)
'''
import argparse
from configparser import ConfigParser
from datetime import datetime
import errno
import os
from pathlib import Path
import re
from shutil import which
import subprocess

from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

from version import VERSION

# Default parameters
FLANK = 500
FUSION = False
TARGET_LEAD = 50
NR_PRIMERS = 4
ALLOWED_MISMATCHES = 5
MAX_MAPPINGS = 5
REFERENCE = None  # depends on the chosen reference genome
SNP = None  # depends on the chosen reference genome

# font file in static dir
FONT = TTFont(
    'mono', (
        f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}'
        '/static/LiberationMono-Regular.ttf'
    )
)

# get thermoparams dir from primer3 build location
try:
    THERMO_PARAMS = f'{Path(which("primer3_core")).parents[0]}/primer3_config/'
except TypeError:
    # will raise TypeError from Path if not valid
    print('primer3 missing from path.')

TMP_FILES = []

colours = [[255, 0, 0],  # red
           [0, 255, 0],  # green
           [0, 0, 255],  # blue
           [255, 0, 255],  # Pink
           [0, 255, 255],  # cyan
           [255, 255, 0]]  # yellow


class Fusion():
    """
    Fusion primer designing related functions

    - split_input(): creates dicts for given coords, used to add attributes
    - fetch_seqs(): calls Sequence.fetch_region() to get sequence from ref
    - flip_fusion_seq(): returns reverse complement where seq is on same
        strand either side of breakpoint
    """

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

            # replacing the a/b with symbols, cos i wrote the whole
            # script with symbols and was lazy to replace them
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
            - coord_dict (dict): nested dict of coordinate dicts
        Returns:
            - coord_dict (dict): nested dict of coordinate dicts with
                the nucleotide sequence added
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
                # invalid side symbol passed
                raise ValueError((
                    "The fusion sequence was marked incorrectly. \n The "
                    f"signs to be used: <>. The sign used:  {dict['SIDE']}"
                ))

        return coord_dict

    def flip_fusion_seq(self, seqs_dict):
        """
        Loops over lsit of dicts of fusion sequences, flips a sequence
        in case two same sides relative to breakpoint are added.

        Args:
            - seqs_dict (dict): dict of dicts for each sequence
        Returns:
            - target_sequence (str): both target sequences joined
            - marked_sequence (str): both marked sequences joined
            - tagged_sequence (str): both tagged sequences joined
        """
        complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            '<': '>', '>': '<', '[': ']', ']': '['
        }

        # DARK_SIDE is required to store original SIDE value, which is
        # used in pick_best_primer() as well as printing coordinates in
        # the right order

        for i, _ in enumerate(seqs_dict):
            # check each dict and tag appropriately
            if seqs_dict[i]['SIDE'] == ">" and seqs_dict[i]['STRAND'] == "1":
                # add DARK_SIDE tag same as side
                seqs_dict[i]['DARK_SIDE'] = ">"

            elif seqs_dict[i]['SIDE'] == "<" and seqs_dict[i]['STRAND'] == "1":
                # add DARK_SIDE tag same as side
                seqs_dict[i]['DARK_SIDE'] = "<"

            elif seqs_dict[i]['SIDE'] == ">" and seqs_dict[i]['STRAND'] == "-1":
                # on rev strand, get reverse complement of nroaml and marked
                # sequence, reverse tagged seq and add DARK_SIDE
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])

                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])

                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = "<"

            elif seqs_dict[i]['SIDE'] == "<" and seqs_dict[i]['STRAND'] == "-1":
                # on rev strand, get reverse complement of nroaml and marked
                # sequence, reverse tagged seq and add DARK_SIDE
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])

                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])

                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = ">"

        if seqs_dict[0]['DARK_SIDE'] == ">" and seqs_dict[1]['DARK_SIDE'] == "<":
            # one set seqs before and after breakpoint => valid (?)
            target_sequence = ''.join([seqs_dict[0]['SEQ'], seqs_dict[1]['SEQ']])
            marked_sequence = ''.join([seqs_dict[0]['MSEQ'], seqs_dict[1]['MSEQ']])
            tagged_string = ''.join([seqs_dict[0]['TSEQ'], seqs_dict[1]['TSEQ']])

        elif seqs_dict[0]['DARK_SIDE'] == "<" and seqs_dict[1]['DARK_SIDE'] == ">":
            # one set seqs before and after breakpoint => valid (?)
            target_sequence = ''.join([seqs_dict[1]['SEQ'], seqs_dict[0]['SEQ']])
            marked_sequence = ''.join([seqs_dict[1]['MSEQ'], seqs_dict[0]['MSEQ']])
            tagged_string = ''.join([seqs_dict[1]['TSEQ'], seqs_dict[0]['TSEQ']])

        else:
            # I think it gets here if both of the coordinates are either
            # before or after the breakpoint?
            raise ValueError((
                "An error occured in designing fusion primers, this design "
                "does not seem possible as both sequences are on the same "
                "side of the breakpoint. Please review the parameters passed."
            ))

        return target_sequence, tagged_string, marked_sequence


class Sequence():
    """
    Sequence related functions

    - fetch_region(): uses samtools to get nucleotide seq. from given reference
    - markup_sequence(): annotates seq., includes marking repeats and snps
    - markup_repeats(): marks repeats longer than 6 bases
    - markup_SNPs(): marks known SNPs in given dbSNP file
    - fetch_known_SNPs(): retrieves kown SNPs from dbSNP file for given region
    """

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
        # call samtools to get sequence from fasta
        cmd = f"samtools faidx {REFERENCE}  {chrom}:{start}-{end}"

        output = subprocess.run(
            cmd, shell=True, check=True, stdout=subprocess.PIPE
        )

        stdout = output.stdout
        stderr = output.stderr

        if stderr:
            stderr = stderr.decode()
            raise RuntimeError(
                "Error in running samtools faidx. Error:\n"
                f"{stderr}"
            )
        else:
            stdout = stdout.decode()

        sequence = ''.join([
            line for line in (stdout.split("\n")) if not re.match('>', line)
        ])

        return sequence

    def markup_sequence(
        self, flank, sequence, FUSION=False, chrom=None, startpos=None,
        endpos=None
    ):
        """
        Marks the target sequence with "-*-", creates a string of
        characters representing the SNPs, target lead (50 bases up and
        down the position)

        Args:
            - flank (int): flank to add around given position
            - sequence (list): list of nucleotide sequences (?)
            - FUSION (bool): bool fusion flag
            - chrom (str): chrom of given sequence
            - startpos (int): start of given sequence
            - endpos (int): end of given sequence
        Returns:
            - sequence_list (list): list of nucleotide sequences (?)
            - tagged_string (str): sequence with snp tags
        """
        if FUSION:
            # loop in main() passes single dicts for each sequence here
            # in case of fusions
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

                sequence_list[TARGET_LEAD - 1] = (
                    f'{sequence_list[TARGET_LEAD - 1]} ] '
                )

            if side == ">":
                startpos = int(sequence['POS']) - FLANK
                endpos = int(sequence['POS'])

                tags[-1] = '*'

                for x in range(0, len(tags)):
                    if x in range(len(tags) - TARGET_LEAD, len(tags) - 1):
                        tags[x] = '-'

                sequence_list[len(tags) - TARGET_LEAD] = (
                    f'{sequence_list[len(tags) - TARGET_LEAD]} ['
                )

            dbSNPs = self.fetch_known_SNPs(SNP, chrom, startpos, endpos)
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
                if (
                    x in range(start - TARGET_LEAD, start)
                    or x in range(end, end + TARGET_LEAD)
                ):
                    tags[x] = '-'

            # Tag the target sequence (50 bases up and down the central
            # nucleotide), regardless of the nature of the variant only
            # one (1) base is tagged as the target.
            sequence_list[flank - TARGET_LEAD] = ' [' + \
                sequence[flank - TARGET_LEAD]

            sequence_list[(- flank + TARGET_LEAD) - 1] = sequence_list[
                (- flank + TARGET_LEAD) - 1] + '] '

            dbSNPs = self.fetch_known_SNPs(
                SNP, chrom, startpos - flank, endpos + flank
            )

            tags = self.markup_repeats(tags, sequence)

            sequence_list, tagged_string = self.markup_SNPs(
                dbSNPs, sequence_list, tags, startpos, endpos
            )

        return sequence_list, tagged_string

    def markup_repeats(self, tags, sequence):
        """
        Marks the repeats in the sequence longer than 6 bases.

        Args:
            - tags (list): empty list to add tags to (?)
            - sequence (str): nucleotide sequence to identify repeats on
        Returns:
            - tags (list): list of added tags
        """
        regex = r"(.{1,})\1\1\1\1\1+"
        matches = re.finditer(regex, sequence, re.MULTILINE)

        for matchNum, match in enumerate(matches, start=1):
            print(
                "Match {matchNum} was found at {start}-{end}: {match}".format(
                    matchNum=matchNum, start=match.start(),
                    end=match.end(), match=match.group()
                )
            )

            for x in range(match.start(), match.end()):
                if len(match.group()) >= 10:
                    tags[x] = '^'
                elif len(match.group()) < 10:
                    tags[x] = '+'

        return tags

    def markup_SNPs(
        self, dbSNPs, sequence_list, tags, startpos, endpos, FUSION=False,
        side=None
    ):
        """
        Marks the known snps depending on the variables passed

        Args:
            - dbSNPs (list): list of snps in given region
            - sequence_list (list): sequnce to annotate with snps
            - tags
            - startpos (int):
            - endpos (int):
            - FUSION (bool):
            - side (str): '>' / '<', indicates side of breakpoint for fusions
        Returns:
            - sequence_list (str): str of joined nucleotide sequences
            - tagged_string (str): sequence but with snp tags
        """

        for dbSNP in dbSNPs:
            # snp records are list of lists, loop over and parse each

            if len(dbSNP) < 8:
                # some field missing or empty line etc.
                continue

            # split out required fields from snp record in vcf
            _, snp_pos, _, snp_ref, _, _, _, _ = dbSNP
            snp_pos = int(snp_pos)

            if FUSION:
                if side == "<":
                    if snp_pos >= startpos + TARGET_LEAD:
                        mask_pos = snp_pos - startpos
                    else:
                        continue
                elif side == ">":
                    if snp_pos <= endpos - TARGET_LEAD:
                        mask_pos = snp_pos - startpos
                    else:
                        continue
            else:
                if (
                    snp_pos >= startpos - TARGET_LEAD
                    and snp_pos <= endpos + TARGET_LEAD
                ):
                    continue

                mask_pos = snp_pos - (startpos - FLANK)

            for i in range(0, len(snp_ref)):
                # should normally just be one base but can be larger

                if len(sequence_list) <= mask_pos + i:
                    # position outside of sequence region
                    break

                # If already masked skip masking it again.
                if (
                    re.search('<', sequence_list[mask_pos + i])
                    or re.search(r'\[', sequence_list[mask_pos + i])
                ):
                    continue

                # add <> around marked bases in sequence, used to make sure
                # it doesn't get masked again
                sequence_list[mask_pos + i] = (
                    f' <{sequence_list[mask_pos + i]}> '
                )
                # tag for colouring snps
                tags[mask_pos + i] = 'X'
            else:
                pass

        sequence_list = "".join(sequence_list).replace(' ', '')
        tagged_string = "".join(tags)

        return sequence_list, tagged_string

    def fetch_known_SNPs(self, tabix_file, chrom, start, end):
        """
        Retrieve SNPs from dbSNP file for given region

        Args:
            - tabix_file (file): file of snps
            - chrom (str): chrom of query region
            - start (int): start coord of query region
            - end (int): end coord of query region
        Returns:
            - snps (list): list of snps within given region
        """
        # call tabix from samtools to get snps in given region from vcf
        cmd = f"tabix {tabix_file}  {chrom}:{start}-{end}"

        output = subprocess.run(
            cmd, shell=True, check=True, stdout=subprocess.PIPE
        )

        stdout = output.stdout
        stderr = output.stderr

        if stderr:
            stderr = stderr.decode()
            raise RuntimeError(
                "Error in retrieving SNPs. Error:\n"
                f"{stderr}"
            )
        else:
            stdout = stdout.decode()

        # build list of lists of snps from output
        snps = [
            line.split("\t") for line in stdout.split("\n")
        ]

        return snps


class Primer3():
    """
    Functions for calling & handling output from primer3

    - template(): just a big old ugly string with config for primer3
    - write_primer3_file(): creates req. .primer3 file to run primer3 with
    - run_primer3(): calls primer3 and generates dict of sequences
    - check_primers(): generates fasta file of primer seqs. and calls smalt
        to map these to given ref., returns dict of seqs. & mapping summary
    - merge_dict(): flattens a nested dict
    - primer_report(): adds summary report for each primer in dict
    - extract_passed_primer_seqs(): returns list with names and seq. of primers
    """

    def template(self, seq_id, seq, flank, len):
        """
        Returns str of params for .primer3 file to run primer3.
        In separate function as it is long and ugly.

        Args:
            - seq_id (str):
            - seq (str):
            - flank (int):
            - len (int):
        Returns:
            - template (str): formatted str of template parameters
        """
        template = f'''SEQUENCE_ID={seq_id}
            SEQUENCE_TEMPLATE={seq}
            SEQUENCE_TARGET={flank},{len}
            PRIMER_FIRST_BASE_INDEX=1
            PRIMER_TASK=generic
            PRIMER_MIN_THREE_PRIME_DISTANCE=3
            PRIMER_EXPLAIN_FLAG=1
            PRIMER_MAX_LIBRARY_MISPRIMING=12.00
            PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
            PRIMER_PRODUCT_SIZE_RANGE=400-800
            PRIMER_NUM_RETURN=5
            PRIMER_MAX_END_STABILITY=9.0
            PRIMER_MAX_SELF_ANY_TH=45.00
            PRIMER_MAX_SELF_END_TH=35.00
            PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00
            PRIMER_PAIR_MAX_COMPL_END_TH=35.00
            PRIMER_MAX_HAIRPIN_TH=24.00
            PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
            PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
            PRIMER_MIN_SIZE=18
            PRIMER_OPT_SIZE=20
            PRIMER_MAX_SIZE=25
            PRIMER_MIN_TM=58.0
            PRIMER_OPT_TM=60.0
            PRIMER_MAX_TM=62.0
            PRIMER_PAIR_MAX_DIFF_TM=5.0
            PRIMER_TM_FORMULA=1
            PRIMER_SALT_MONOVALENT=50.0
            PRIMER_SALT_CORRECTIONS=1
            PRIMER_SALT_DIVALENT=1.5
            PRIMER_DNTP_CONC=0.6
            PRIMER_DNA_CONC=50.0
            PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
            PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
            PRIMER_LOWERCASE_MASKING=0
            PRIMER_MIN_GC=30.0
            PRIMER_MAX_GC=70.0
            PRIMER_MAX_NS_ACCEPTED=0
            PRIMER_MAX_POLY_X=4
            PRIMER_OUTSIDE_PENALTY=0
            PRIMER_GC_CLAMP=0
            PRIMER_LIBERAL_BASE=1
            PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
            PRIMER_PICK_ANYWAY=1
            PRIMER_WT_TM_LT=1.0
            PRIMER_WT_TM_GT=1.0
            PRIMER_WT_SIZE_LT=1.0
            PRIMER_WT_SIZE_GT=1.0
            PRIMER_WT_GC_PERCENT_LT=0.0
            PRIMER_WT_GC_PERCENT_GT=0.0
            PRIMER_WT_SELF_ANY_TH=0.0
            PRIMER_WT_SELF_END_TH=0.0
            PRIMER_WT_HAIRPIN_TH=0.0
            PRIMER_WT_NUM_NS=0.0
            PRIMER_WT_LIBRARY_MISPRIMING=0.0
            PRIMER_WT_SEQ_QUAL=0.0
            PRIMER_WT_END_QUAL=0.0
            PRIMER_WT_POS_PENALTY=0.0
            PRIMER_WT_END_STABILITY=0.0
            PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
            PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
            PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
            PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
            PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
            PRIMER_PAIR_WT_DIFF_TM=0.0
            PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
            PRIMER_PAIR_WT_COMPL_END_TH=0.0
            PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
            PRIMER_PAIR_WT_PR_PENALTY=1.0
            PRIMER_PAIR_WT_IO_PENALTY=0.0
            PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
            PRIMER_INTERNAL_WT_SIZE_LT=1.0
            PRIMER_INTERNAL_WT_END_QUAL=0.0
            PRIMER_INTERNAL_MAX_SELF_END=12.00
            PRIMER_QUALITY_RANGE_MIN=0
            PRIMER_PAIR_MAX_COMPL_END=3.00
            PRIMER_PRODUCT_MAX_TM=1000000.0
            PRIMER_INTERNAL_MAX_SIZE=27
            PRIMER_INTERNAL_WT_SELF_ANY=0.0
            PRIMER_INTERNAL_MAX_POLY_X=5
            PRIMER_INTERNAL_WT_SIZE_GT=1.0
            PRIMER_SEQUENCING_ACCURACY=20
            PRIMER_INTERNAL_WT_TM_GT=1.0
            PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
            PRIMER_INTERNAL_MAX_GC=80.0
            PRIMER_PAIR_WT_COMPL_ANY=0.0
            PRIMER_PICK_INTERNAL_OLIGO=0
            PRIMER_MAX_SELF_END=3.00
            PRIMER_QUALITY_RANGE_MAX=100
            PRIMER_INTERNAL_DNTP_CONC=0.0
            PRIMER_INTERNAL_MIN_SIZE=18
            PRIMER_INTERNAL_MIN_QUALITY=0
            PRIMER_SEQUENCING_INTERVAL=250
            PRIMER_INTERNAL_SALT_DIVALENT=1.5
            PRIMER_MAX_SELF_ANY=8.00
            PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
            PRIMER_PAIR_WT_COMPL_END=0.0
            PRIMER_INTERNAL_OPT_TM=60.0
            PRIMER_SEQUENCING_SPACING=500
            PRIMER_INTERNAL_MAX_SELF_ANY=12.00
            PRIMER_MIN_END_QUALITY=0
            PRIMER_INTERNAL_MIN_TM=57.0
            PRIMER_PAIR_MAX_COMPL_ANY=8.00
            PRIMER_SEQUENCING_LEAD=50
            PRIMER_PICK_LEFT_PRIMER=1
            PRIMER_INTERNAL_OPT_SIZE=20
            PRIMER_WT_TEMPLATE_MISPRIMING=0.0
            PRIMER_MAX_END_GC=5
            PRIMER_MIN_QUALITY=0
            PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
            PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
            PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
            PRIMER_WT_SELF_ANY=0.0
            PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
            PRIMER_INTERNAL_WT_NUM_NS=0.0
            PRIMER_INTERNAL_WT_SELF_END=0.0
            PRIMER_PRODUCT_OPT_SIZE=0
            PRIMER_PRODUCT_OPT_TM=0.0
            PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
            PRIMER_INSIDE_PENALTY=-1.0
            PRIMER_INTERNAL_MIN_GC=20.0
            PRIMER_PRODUCT_MIN_TM=-1000000.0
            PRIMER_INTERNAL_SALT_MONOVALENT=50.0
            PRIMER_WT_SELF_END=0.0
            PRIMER_INTERNAL_DNA_CONC=50.0
            PRIMER_PICK_RIGHT_PRIMER=1
            PRIMER_INTERNAL_MAX_TM=63.0
            PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
            PRIMER_INTERNAL_WT_TM_LT=1.0
            PRIMER_THERMODYNAMIC_PARAMETERS_PATH={THERMO_PARAMS}
            ='''

        # stupid regex to remove leading tabs but it works so...
        template = re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', template, flags=re.M)

        return template

    def write_primer3_file(self, seq_id, seq, primer3_file):
        """
        Creates a *.primer3 file required for PRIMER3 to generate primers
        for a given sequence

        Args:
            - seq_id (int):
            - seq (str):
            - primer3_file (str): file name for primer3 file

        Returns: None

        Outputs: primer3 run file
        """
        upstream_end = 1
        downstream_start = 2

        for i in range(0, len(seq)):
            if seq[i] == '[':
                upstream_end = i
            elif seq[i] == ']':
                downstream_start = i

        flank = upstream_end
        length = downstream_start - upstream_end + 1
        template = self.template(seq_id, seq, flank, length)

        with open(primer3_file, 'w+') as outfile:
            outfile.write(template)

        outfile.close()

    def run_primer3(self, seq_id, seq, primer3_file=""):
        """
        Creates a primer3 file. Runs PRIMER3 and captures its output and
        returns it as a dictionary where the name of the primer is the
        key and the sequence is the value

        Args:
            -
        Returns:
            -
        """
        if primer3_file == "":
            primer3_file = seq_id + ".primer3"

        # prepend absolute path to file
        primer3_file = f"{Path(__file__).parent.absolute()}/{primer3_file}"

        # generate required .primer3 file with config in tmp dir
        primer3_file = re.sub("[<>:]", "_", primer3_file)
        TMP_FILES.append(primer3_file)
        self.write_primer3_file(seq_id, seq, primer3_file)

        # call primer3 to generate primers
        cmd = f"primer3_core {primer3_file}"

        output = subprocess.run(
            cmd, shell=True, check=True, stdout=subprocess.PIPE
        )

        stdout = output.stdout
        stderr = output.stderr

        if stderr:
            stderr = stderr.decode()
            raise RuntimeError(
                "Error in running primer3 to design primers. Error:\n"
                f"{stderr}"
            )
        else:
            stdout = stdout.decode()

        output_dict = dict()

        for line in stdout.split('\n'):
            if line == '=':
                break

            key, value = line.split("=")
            output_dict[key] = value

        return output_dict

    def check_primers(self, region_id, target_region, primer3_dict,
                      chromo, startpos, endpos, FUSION=False, seqs=None):
        """
        Creates a fasta file with sequences of the primers, which are
        requrired by SMALT to map to the reference genome. The smalt result
        is saved in a nested dictionary containing mapping summary of each
        primer.

        Args:
            - region_id (str): name of region used for naming smalt files
            - target_region (str): nucleotide sequence of ref
            - primer3_dict (dict): output from primer3 with designed primers
            - chrom (str): chromsome of primer
            - startpos (int): start position to write in report
            - endpos (int): end position to write in report
            - FUSION (bool): if design is a Fusion
            - seqs (str): dict of nucleotide sequences and positions
        Returns:
            - seq_dict (dict): dictionary of sequences, positions and smalt
                results (total mappings etc.)
        """
        region_id = re.sub("[<>:]", "_", region_id)

        primers_file = region_id + "_p3seq.fasta"
        TMP_FILES.append(primers_file)

        with open(primers_file, 'w+') as primerfasta:
            primerfasta.write(">FULLSEQ\n" + target_region + "\n")

            for i in range(0, 5):
                pl_id = "PRIMER_LEFT_{}_SEQUENCE".format(i)
                pr_id = "PRIMER_RIGHT_{}_SEQUENCE".format(i)

                if pr_id not in primer3_dict or pl_id not in primer3_dict:
                    continue

                primerfasta.write(f">{pl_id}\n{primer3_dict[pl_id]}\n")
                primerfasta.write(f">{pr_id}\n{primer3_dict[pr_id]}\n")

            primerfasta.close()

        smalt_out = region_id + ".smalt"
        TMP_FILES.append(smalt_out)

        # .sma index for reference generated by smalt sometimes has .fasta
        # extension removed, check for both and raise error if can't find one
        ref_file = os.path.expanduser(REFERENCE).rstrip(".fasta")

        if Path(f'{ref_file}.sma').is_file():
            # index with no .fasta suffix
            pass
        elif Path(f'{ref_file}.fasta.sma').is_file():
            # index with .fasta suffix
            ref_file = f'{ref_file}.fasta'
        else:
            print((
                'No .sma file found for reference, has it first been indexed '
                'by smalt?'
            ))
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                f'{ref_file}.fasta OR {ref_file}.fasta.sma'
            )

        # run smalt to map primer sequences to reference
        cmd = f"smalt map -x -d -1 {ref_file} {primers_file} > {smalt_out}"
        subprocess.call(cmd, shell=True)

        seq_dict = {}

        with open(smalt_out, "r") as f:
            for line in f.readlines():
                if line.startswith("@") or line.startswith("FULLSEQ"):
                    continue

                # split line & parse
                line = line.split('\t')
                name = line[0]
                chrom = line[2]
                pos = line[3]
                mismatch = len(line[9]) - int(line[12].split(':')[2])

                if mismatch > ALLOWED_MISMATCHES:
                    # too many mismatches, drop
                    print(name)
                    print(mismatch)
                    print('')
                    continue

                if name not in seq_dict:
                    seq_dict[name] = {}
                    seq_dict[name]['CHR'] = [chrom]
                    seq_dict[name]['POS'] = [pos]
                    seq_dict[name]['MISMATCH'] = [str(mismatch)]

                else:
                    seq_dict[name]['CHR'].append(chrom)
                    seq_dict[name]['POS'].append(pos)
                    seq_dict[name]['MISMATCH'].append(str(mismatch))

        left_dict = {}
        right_dict = {}

        if FUSION:
            # as there might be 2 different sequences on 2 chromosomes the
            # dictionary has to be broken down and summaries have to be made
            # separately of each other

            for key in seq_dict:
                if "RIGHT" in key:
                    right_dict[key] = seq_dict[key]
                elif "LEFT" in key:
                    left_dict[key] = seq_dict[key]

            for _, dicts in seqs.items():
                if dicts["DARK_SIDE"] == ">":
                    chromo = dicts['CHR']
                    startpos = int(dicts['POS']) - FLANK
                    endpos = int(dicts['POS'])
                    left_dict = self.primer_report(
                        chromo, startpos, endpos, left_dict
                    )

                elif dicts["DARK_SIDE"] == "<":
                    chromo = dicts['CHR']
                    startpos = int(dicts['POS'])
                    endpos = int(dicts['POS']) + FLANK
                    right_dict = self.primer_report(
                        chromo, startpos, endpos, right_dict
                    )

            # merge two dictionaries back into one
            seq_dict = self.merge_dicts(left_dict, right_dict)
        else:
            seq_dict = self.primer_report(chromo, startpos, endpos, seq_dict)

        # check if seq_dict is empty and raise error
        if not seq_dict:
            raise RuntimeError(
                "Error in generating primers - no primers generated by primer3"
                " OR no mappings output from smalt. Exiting now."
            )

        return seq_dict

    def merge_dicts(self, *dicts):
        """
        Given any number of dicts, shallow copy and merge into a new dict,
        precedence goes to key value pairs in latter dicts.

        Args:
            - *dict_args(dict): multiple dicts
        Returns:
            - merged_dicts (dict): dict of merged dicts
        """
        merged_dicts = {}
        for dictionary in dicts:
            merged_dicts.update(dictionary)

        return merged_dicts

    def primer_report(self, chromo, startpos, endpos, seq_dict):
        """
        Creating a nested dictionary with summary report for each primer

        Args:
            - chrom (str): chromosome of primer
            - startpos (int): start position of primer mapping
            - endpos (int): end position of primer mapping
            - seq_dict (dict): dict of nucleotide sequence, coordinates and
                mapping info etc.
        Returns:
            - seq_dict (dict): dict of nucleotide sequence, coordinates,
                mapping info and added summary text for displaying in report
        """
        for k, v in seq_dict.items():
            # dict of dicts with primer information, each sub dict has
            # chr, pos and list of mismatches
            uniq_chr = []
            mis_chr = []

            for pos, chrom in enumerate(seq_dict[k]['CHR']):
                # same chromosome, position in range and no mismatch
                if (
                    seq_dict[k]['CHR'][pos] == chromo
                    and int(seq_dict[k]['POS'][pos])
                    in range(startpos - FLANK, endpos + 1 + FLANK)
                    and seq_dict[k]['MISMATCH'][pos] == '0'
                ):
                    # uniquely mapped primer and in within region
                    uniq_chr.append(seq_dict[k]['CHR'][pos])
                else:
                    # not uniquely mapped / outside region
                    mis_chr.append(seq_dict[k]['CHR'][pos])

            if len(uniq_chr) == 1 and len(mis_chr) == 0:
                # uniquely mapped
                seq_dict[k]['MAPPING_SUMMARY'] = 'unique mapping'

            elif len(uniq_chr) == 0 and len(mis_chr) == 1:
                # unique but chromosome doesnt match, somehow primer has
                # been designer but seems to not map back to target region
                seq_dict[k]['MAPPING_SUMMARY'] = 'unique - wrong chromosome'

            elif len(uniq_chr) + len(mis_chr) in range(2, 6):
                # 2-5 mismatches
                seq_dict[k]['MAPPING_SUMMARY'] = (
                    f'{len(uniq_chr + mis_chr)} mapping(s) on chr '
                    f'{", ".join(uniq_chr + mis_chr)} with '
                    f'{", ".join(seq_dict[k]["MISMATCH"])} mismatches'
                )
            elif len(uniq_chr) + len(mis_chr) > 5:
                # many mismatches
                seq_dict[k]['MAPPING_SUMMARY'] = (
                    f'{len(uniq_chr + mis_chr)} mapping(s)'
                )
            else:
                # shouldn't get here, added an else to stop a key error later
                seq_dict[k]['MAPPING_SUMMARY'] = (
                    'Error in primer, could not determine total mismatches'
                )

        return seq_dict

    def extract_passed_primer_seqs(self, primer3_results, passed_primers):
        """
        Extracting the name and sequnece of the primers generated by primer3

        Args:
            - primer3_results (dict): output from primer3
            - passed_primers (dict): dictionary of sequences, positions and
                smalt results (total mappings etc.)
        Returns:
            - primer_seqs (list): list of lists with primer name and sequence
        """
        primer_seqs = []

        for primer in passed_primers:
            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            primer_seqs.append([name, primer3_results[primer]])

        return primer_seqs


class Report():
    """
    Functions to generate PDF report. Some really weird f strings for
    adding appropriate padding for displaying report.

    - make_primer_mapped_str(): formats sequence & primers for displaying
        in report
    - align_primers_to_seq(): checks if any primers match any regions of the
        sequence of interest?
    - check_if_primer_clash(): check if any primers clash with bases between
        start and end coords
    - revDNA(): returns reverse compliment of given nucleotide sequence

    - pretty_pdf_primer_data(): formats header info for pdf
    - pretty_pdf_fusion_mappings(): annotates the nested dict of sequences
        for adding to the report
    - pretty_pdf_mappings(): adds the mapped sequence to the pdf report
    - pretty_pdf_method(): writes the blurb from method_blurb() to the report
    - method_blurb(): generates blurb text for report, inc what files used etc.
    - pretty_primer_data(): generates the .txt file if specified in args
    """

    def make_primer_mapped_str(self, target_sequence, passed_primer_seqs):
        """
        Creates a mapping string out of the target sequence and primers.
        Used to be printed on the PDF report.

        Args:
            - target_sequence (str):
            - passed_primer_seqs (dict):
        Returns:
            - mapped_strings (list):
            - mapped_colours (list):
        """
        mappings = self.align_primers_to_seq(
            target_sequence, passed_primer_seqs
        )

        mapped_strings = []
        mapped_strings.append([" "] * len(target_sequence))

        mapped_colours = []
        mapped_colours.append([-1] * len(target_sequence))

        for mapping in mappings:
            name, primer, pos, strand = mapping

            primer_nr = int(re.sub(r'.*_(\d)', r'\1', name))
            mapping_index = 0

            arrows = int((len(primer) - len(name) - 2) / 2)
            arrow_type = ">"

            if strand == 1:
                # if the primer is on the reverse strand
                primer = self.revDNA(primer)
                arrow_type = "<"

            tag = arrow_type * arrows + " " + name + " " + arrow_type * arrows

            if len(tag) < len(primer):
                tag += arrow_type

            primer = list(tag)
            scan_index = 0

            while(1):
                if self.check_if_primer_clash(
                        mapped_strings[scan_index],
                        pos,
                        pos + len(primer)):
                    scan_index += 1

                    if len(mapped_strings) >= scan_index:
                        mapped_strings.append([" "] * len(target_sequence))
                        mapped_colours.append([-1] * len(target_sequence))
                else:
                    break

            mapping_index = scan_index

            for i in range(0, len(primer)):
                mapped_strings[mapping_index][pos + i] = primer[i]
                mapped_colours[mapping_index][pos + i] = primer_nr

        for i in range(0, len(mapped_strings)):
            mapped_strings[i] = "".join(mapped_strings[i])

        return mapped_strings, mapped_colours

    def check_if_primer_clash(self, mapped_list, start, end):
        """
        Checks whether the primers clash

        Args:
            - mapped_list (list): list of mapped positions
            - start (int): start pos
            - end (int): end pos
        Returns:
            - bool: True if present in list else False
        """
        for i in range(start, end):
            if mapped_list[i] != " ":
                return True

        return False

    def align_primers_to_seq(self, seq, all_primers):
        """
        Checking if the primers and their reverse sequences match any
        regions of the sequence of interest adding everything to a list.
        Used in make_primers_mapped_strings().

        Args:
            - seq ():
            - all_primers():
        Returns:
            - mappings ():
        """
        primers = []
        rev_primers = []

        for primer_set in all_primers:
            name, primer = primer_set
            primers.append([name, primer])
            rev_primers.append([name, self.revDNA(primer)])

        mappings = []

        for i in range(0, len(seq)):
            # for each primer checks at each base of target sequence if the
            # following bases match the primer sequence, super crude way of
            # 'mapping' the primer seuqnece to know where in the target
            # region it is for adding <<< to the report
            for primer_set in primers:
                (name, primer) = primer_set
                primer_len = len(primer)

                if "".join(seq[i:i + primer_len]) == primer:
                    mappings.append([name, primer, i, 0])

            for primer_set in rev_primers:
                name, primer = primer_set
                primer_len = len(primer)

                if "".join(seq[i:i + primer_len]) == primer:
                    mappings.append([name, primer, i, 1])

        return mappings

    def revDNA(self, string):
        """
        Returns reverse compliment of given nucleotide sequence.

        Args:
            - string (str): nucleotide sequence
        Returns:
            - rev_str (str): nucleotide sequence but magically reversed
        """
        rev_bases = {
            'A': 'T', 'a': 'T', 'C': 'G', 'c': 'G',
            'G': 'C', 'g': 'C', 'T': 'A', 'T': 'A',
            '-': '-'
        }

        rev_str = len(string) * [None]

        for i in range(0, len(string)):
            rev_str[len(string) - i - 1] = rev_bases[string[i]]

        rev_str = "".join(rev_str)

        return rev_str

    def pretty_pdf_primer_data(
            self, c, y_offset, primer3_results, passed_primers, width, args,
            FUSION=False, chrom=None, startpos=None, endpos=None,
            coord_dict=None):
        """
        Adds the header primer information to the output PDF file

        Args:
            -
        Returns:
            -
        """
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8

        if args.grch37:
            ref = 'GRCh37'
        else:
            ref = 'GRCh38'

        if FUSION:
            # fusion header
            c.drawString(
                40, y_offset,
                (
                    f"Primer design report for a fusion between chr: "
                    f"{coord_dict[0]['CHR']} position: {coord_dict[0]['POS']} "
                    f"and chr: {coord_dict[1]['CHR']} position: "
                    f"{coord_dict[1]['POS']} - ({ref})"
                )
            )
        else:
            # normal primer design header
            if startpos == endpos:
                c.drawString(
                    40,
                    y_offset, (
                        f"Primer design report for chr {chrom} position: "
                        f"{startpos} - ({ref})"
                    )
                )
            else:
                c.drawString(
                    40, y_offset, (
                        f"Primer design report for chr: {chrom} range: "
                        f"{startpos}-{endpos} - ({ref})"
                    )
                )

        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 16

        # weird spacing to align in report
        c.drawString(
            40, y_offset,
            (
                f"ID{' ' * 9}%GC{' ' * 4}TM{' ' * 5}Primer sequence{' ' * 23}"
                "Mapping(s)"
            )
        )

        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8

        for primer in sorted(passed_primers):
            # loop over dict of primers
            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            if name == "RIGHT_0":
                y_offset -= 8

            # set vars to write into report
            gc = float(primer3_results['PRIMER_' + name + '_GC_PERCENT'])
            tm = float(primer3_results['PRIMER_' + name + '_TM'])
            seq = primer3_results['PRIMER_' + name + '_SEQUENCE']
            summary = passed_primers[
                'PRIMER_' + name + '_SEQUENCE']['MAPPING_SUMMARY']

            # more weird spacing for report alignment, adds in primer details
            c.drawString(
                40, y_offset,
                (
                    f"{'':10} {gc:.2f}  {tm:.2f}  "
                    f"{seq:25}{'':13} {summary}{'':12}"
                )
            )

            primer_nr = re.sub(r'.*_(\d)', r'\1', name)

            c.setFillColorRGB(
                colours[int(primer_nr)][0], colours[int(primer_nr)][1],
                colours[int(primer_nr)][2]
            )

            c.drawString(40, y_offset, name)
            c.setFillColorRGB(0, 0, 0)
            y_offset -= 8

        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8
        y_offset -= 8

        return y_offset

    def pretty_pdf_fusion_mappings(
        self, top_offset, c, coord_dict, passed_primer_seqs, FUSION=False
    ):
        """
        Function to interogate the nested dictionary.
        The first sequence with the position of interest at the end is
        printed. It is followed by the other one.
        DARK_SIDE varible allows to track the change of the side of the
        flipped sequence.
        If sequence of interest was after given position and then
        reverse-complemented the signs would be < > for side and
        dark_side respectively.

        Args:

        Returns:

        """
        # spaces required to match the positions of two breakpoints when
        # pdf is created
        spaces = ' ' * ((FLANK + 1) % 80)

        if FUSION:
            if coord_dict[0]['TSEQ'][-1] == "*":
                target_sequence = coord_dict[0]['SEQ']
                tagged_string = coord_dict[0]['TSEQ']

                if coord_dict[0]['SIDE'] == "<" and coord_dict[0]['STRAND'] == "-1":
                    # the coordinates are printed FLANK bases down the given
                    # position s
                    base1 = int(coord_dict[0]['POS']) + FLANK

                elif coord_dict[0]['SIDE'] == ">" and coord_dict[0]['STRAND'] == "1":
                    base1 = int(coord_dict[0]['POS']) - FLANK

                side = coord_dict[0]['DARK_SIDE']
                darkside = coord_dict[0]['DARK_SIDE']

                primer_strings, primer_colours = self.make_primer_mapped_str(
                    target_sequence, passed_primer_seqs)

                top_offset = self.pretty_pdf_mappings(
                    top_offset, target_sequence, tagged_string, primer_strings,
                    primer_colours, base1, c, side, darkside
                )

                target_sequence = spaces + coord_dict[1]['SEQ']
                tagged_string = spaces + coord_dict[1]['TSEQ']
                base1 = int(coord_dict[1]['POS'])
                side = coord_dict[1]['DARK_SIDE']
                darkside = coord_dict[1]['DARK_SIDE']

                primer_strings, primer_colours = self.make_primer_mapped_str(
                    target_sequence, passed_primer_seqs
                )

                top_offset = self.pretty_pdf_mappings(
                    top_offset, target_sequence, tagged_string, primer_strings,
                    primer_colours, base1, c, side, darkside
                )

            elif coord_dict[1]['TSEQ'][-1] == "*":

                target_sequence = coord_dict[1]['SEQ']
                tagged_string = coord_dict[1]['TSEQ']

                if coord_dict[1]['SIDE'] == "<" and coord_dict[1]['STRAND'] == "-1":
                    # the coordinates are printed FLANK bases down the given
                    # position s
                    base1 = int(coord_dict[1]['POS']) + FLANK

                elif coord_dict[1]['SIDE'] == ">" and coord_dict[1]['STRAND'] == "1":
                    # the coordinates are printed FLANK bases down the given
                    # position
                    base1 = int(coord_dict[1]['POS']) - FLANK

                side = coord_dict[1]['DARK_SIDE']
                darkside = coord_dict[1]['DARK_SIDE']

                primer_strings, primer_colours = self.make_primer_mapped_str(
                    target_sequence, passed_primer_seqs
                )

                top_offset = self.pretty_pdf_mappings(
                    top_offset, target_sequence, tagged_string, primer_strings,
                    primer_colours, base1, c, side, darkside
                )

                target_sequence = spaces + coord_dict[0]['SEQ']
                tagged_string = spaces + coord_dict[0]['TSEQ']
                base1 = int(coord_dict[0]['POS'])
                side = coord_dict[0]['DARK_SIDE']
                darkside = coord_dict[0]['DARK_SIDE']

                primer_strings, primer_colours = self.make_primer_mapped_str(
                    target_sequence, passed_primer_seqs
                )

                top_offset = self.pretty_pdf_mappings(
                    top_offset, target_sequence, tagged_string, primer_strings,
                    primer_colours, base1, c, side, darkside
                )

        return top_offset

    def pretty_pdf_mappings(
            self, top_offset, target_sequence, tagged_string, primer_strings,
            primer_colours, base1, c, side=None, darkside=None):
        """
        Writes the reference sequence and primers with markup to the report

        Args:

        Returns:
        """
        # throughout this mess p_line is the 'primer line' which is actually
        # the reference target sequence region and m_line is the 'markup line'
        # that has the sequence annotation
        spaces = ' ' * ((FLANK + 1) % 80)

        # The flip of one of the sequences defines how the numeration of
        # coordinates go, either in descending or ascending orders

        for i in range(0, len(target_sequence), 80):
            # base1 is initial coordinate, iterating in blocks of 80 bases
            # to the end of the ref sequence. This defines the width to the
            # beginning of the m_line to match the beginning of p_line.
            # All of this basically just aligns both sequence and the markup.
            if FUSION:
                if side == '>' and darkside == '>':
                    # forward strand and before break
                    start_coord = str(base1 + i)
                    m_line_pad = ' ' * len(start_coord)

                    p_line = f"{start_coord}  {target_sequence[i: i + 80]}"
                    m_line = f"{m_line_pad}  {tagged_string[i: i + 80]}"

                elif side == '<' and darkside == '<':
                    # reverse strand and after break
                    start_coord = str(base1 + i - len(spaces))
                    m_line_pad = ' ' * len(start_coord)

                    p_line = f"{start_coord}  {target_sequence[i: i + 80]}"
                    m_line = f"{ m_line_pad}  {tagged_string[i: i + 80]}"

                elif side == '>' and darkside == '<':
                    # forward strand and after break
                    start_coord = str(base1 - i + len(spaces))
                    m_line_pad = ' ' * len(start_coord)

                    p_line = f"{start_coord}  {target_sequence[i: i + 80]}"
                    m_line = f"{m_line_pad}  {tagged_string[i: i + 80]}"

                elif side == '<' and darkside == '>':
                    # reverse strand and before break
                    start_coord = str(base1 - i)
                    m_line_pad = ' ' * len(start_coord)

                    p_line = f"{start_coord}  {target_sequence[i: i + 80]}"
                    m_line = f"{m_line_pad}  {tagged_string[i: i + 80]}"
            else:
                # normal primers
                start_coord = str(base1 + i)
                m_line_pad = ' ' * len(start_coord)

                p_line = f"{start_coord}  {target_sequence[i: i + 80]}"
                m_line = f"{m_line_pad}  {tagged_string[i: i + 80]}"

            x_offset = 40

            for k in range(0, len(p_line)):
                # Setting colour of the SNPs and the indicated position
                if m_line[k] == "X":
                    c.setFillColorRGB(255, 0, 0)
                elif m_line[k] == "*":
                    c.setFillColorRGB(0, 190, 0)
                elif m_line[k] == "+":
                    c.setFillColorRGB(0, 0, 255)
                elif m_line[k] == "^":
                    c.setFillColorRGB(0, 0, 255)

                # write the reference sequence to pdf
                c.drawString(x_offset, top_offset, p_line[k])
                x_offset += stringWidth(" ", 'mono', 8)
                c.setFillColorRGB(0, 0, 0)

            top_offset -= 8

            m_line = re.sub(r'X', r' ', m_line)

            if re.search(r'\*', m_line) or re.search(r'\-', m_line):
                # setting colour for the TARGET LEAD (50 base pairs around the
                # position of interest) & writes target region line to pdf
                c.setFillColorRGB(0, 190, 0)
                c.drawString(40, top_offset, m_line)
                c.setFillColorRGB(0, 0, 0)

                top_offset -= 8

            for j in range(0, len(primer_strings)):
                # write the rest of the markup to the pdf
                primer_string = primer_strings[j]
                primer_colour = primer_colours[j]

                line = primer_string[i: i + 80]
                if re.match(r'^ *$', line):
                    # empty line
                    continue

                # padding width to match sequence line, weird format but works
                x_offset = 40 + stringWidth(" ", 'mono', len(m_line_pad)) * 10

                for k in range(i, i + 80):
                    if k > len(target_sequence) - 1:
                        break

                    if primer_colour[k] >= 0:
                        # set appropriate colouring
                        c.setFillColorRGB(
                            colours[primer_colour[k]][0],
                            colours[primer_colour[k]][1],
                            colours[primer_colour[k]][2]
                        )

                    # write the primer sequence markers (i.e. left & right)
                    c.drawString(x_offset, top_offset, primer_string[k])
                    x_offset += stringWidth(" ", 'mono', 8)
                    c.setFillColorRGB(0, 0, 0)

                top_offset -= 8

            top_offset -= 8

        return top_offset

    def pretty_pdf_method(self, top_offset, args, c, seqs=None):
        """
        Function to do something to the pdf - not sure what...
            - think it might just be writing the blurb to the report?

        Args:

        Returns:
        """
        lines = self.method_blurb(args, seqs)

        top_offset = 140

        for line in lines:
            c.drawString(40, top_offset, line)
            top_offset -= 8

    def method_blurb(self, args, seqs=None):
        """
        Generates summary blurb text to write into report

        Args:
            - args (args): cmd line arguments
            - seqs (dict): dict of sequence info
        Returns:
            - lines (list): lines of summary text
        """

        lines = []

        if FUSION:
            lines.append((
                f"Primer design report for a fusion between "
                f"{seqs[0]['CHR']}:{seqs[0]['POS']} and "
                f"{seqs[1]['CHR']}:{seqs[1]['POS']}"
            ))

            for regionid, region_dict in seqs.items():
                if region_dict['STRAND'] == "-1":
                    reverse_region = region_dict['CHR']
                    lines.append((
                        f'The sequence on chromosome {reverse_region} was '
                        'reverse-complemented to produce this fusion sequence'
                    ))
                else:
                    continue

        if args.grch37:
            ref = 'GRCh37'
        else:
            ref = 'GRCh38'

        # get versions of tools used
        primer3_version = subprocess.run(
            f"primer3_core --about", shell=True, check=True,
            stdout=subprocess.PIPE
        )
        primer3_version = primer3_version.stdout.decode('utf-8').split()[-1]

        samtools_version = subprocess.run(
            f"samtools --version", shell=True, check=True,
            stdout=subprocess.PIPE
        )
        samtools_version = samtools_version.stdout.decode('utf-8').split()[1]

        smalt_version = subprocess.run(
            f"smalt version", shell=True, check=True,
            stdout=subprocess.PIPE
        )
        smalt_version = smalt_version.stdout.decode('utf-8').split()[8]

        # record cmd line args used, split into 3 per line to limit length
        cmd_args = []
        [cmd_args.append(f'--{k} {v}') for k, v in vars(args).items() if v][0]
        cmd_args = ['; '.join(cmd_args[i:i + 3]) for i in range(0, len(cmd_args), 3)]

        lines.append('')
        lines.append(
            f'Created at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}'
        )
        lines.append(f'SNP file used: {Path(SNP).stem}')
        lines.append(f'Human reference version: {ref} ({Path(REFERENCE).stem})')

        # add final footer text, split over lines as line breaks refused
        # to work and I don't care enough to debug reportlab
        lines.append('')
        lines.append((
            'Common SNP annotation: A common SNP is one that has allele frequency (AF) '
            'higher than or equal to 1%'
        ))
        lines.append((
            f'in the gnomad database version {SNP_VERSION}.'
        ))

        lines.append('')
        lines.append(f'primer-designer version: {VERSION}')
        lines.append(f'primer3 version: {primer3_version}')
        lines.append(f'samtools version: {samtools_version}')
        lines.append(f'smalt version: {smalt_version}')
        lines.append(f'cmd line args used: {cmd_args[0]}')
        # write extra cmd line args with indent for nice alignment
        [lines.append(f'{" " * 20}{x}') for x in cmd_args[1:]]

        return lines

    def pretty_primer_data(self, outfile, primer3_results, passed_primers,
                           chrom, startpos, endpos, target_sequence=None,
                           FUSION=False, coord_dict=None):
        """
        Function to output the report in txt format if specified.
        No primers would be displayed on the target sequence

        Args:

        Returns: None

        Outputs: {outfile}.txt - text file with primer info
        """
        fh = open(outfile, 'w')

        lines = []
        if FUSION:
            lines.append(
                f"Primer design report for a fusion between chr: "
                f"{coord_dict[0]['CHR']} position: {coord_dict[0]['POS']} "
                f"and chr: {coord_dict[1]['CHR']} position: "
                f"{coord_dict[1]['POS']} "
            )

        elif startpos == endpos:
            lines.append(
                "Primer design report for chr: {} position: {}".format(
                    chr, startpos))

        else:
            lines.append(
                "Primer design report for chr: {} range: {}-{}".format(
                    chr, startpos, endpos
                )
            )

        lines.append("ID\t%GC\tTM\tPrimer sequence\tBest primer\tMapping(s)")

        for primer in sorted(passed_primers):
            if primer == 'FULLSEQ':
                continue

            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            lines.append("\t".join([
                name, primer3_results["PRIMER_" + name + "_GC_PERCENT"],
                primer3_results["PRIMER_" + name + "_TM"],
                primer3_results["PRIMER_" + name + "_SEQUENCE"],
                passed_primers["PRIMER_" + name + "_SEQUENCE"]['MAPPING_SUMMARY']
            ]))

        lines.append("\n")
        lines.append("Consensus sequence:\n")
        lines. append(target_sequence)
        fh.write("\n".join(lines))
        fh.close()


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
            "Design primers around breakpoint. --b1 and --b2 must be specified "
            "and in this format chr:pos:side:strand, where side = a or b "
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
        else:
            REFERENCE = config['REFERENCE']['REF_38']
            SNP = config['REFERENCE']['SNP_37']
    else:
        # no config file, try read from env
        if args.grch37:
            REFERENCE = os.environ.get('REF_37', None)
            SNP = os.environ.get('SNP_37', None)
        else:
            REFERENCE = os.environ.get('REF_38', None)
            SNP = os.environ.get('SNP_37', None)
        VERSION = os.environ.get('PRIMER_VERSION', None)
        SNP_VERSION = os.environ.get('SNP_VERSION', None)

    if not REFERENCE and not SNP:
        raise ValueError('Missing path for reference or snp')

    return REFERENCE, SNP, VERSION, SNP_VERSION


def main():
    global FUSION
    global TARGET_LEAD
    global SNP
    global FLANK
    global NR_PRIMERS
    global ALLOWED_MISMATCHES
    global MAX_MAPPINGS
    global REFERENCE
    global VERSION
    global FONT
    global TMP_FILES
    global SNP_VERSION

    # parse args, load in file paths from config
    args = parse_args()
    REFERENCE, SNP, VERSION, SNP_VERSION = load_config(args)

    # check required tools installed and on PATH
    for tool in ['samtools', 'tabix', 'primer3_core', 'smalt']:
        assert which(tool), f'{tool} is not on path, is it installed?'

    # initialise classes
    fusion = Fusion()
    sequence = Sequence()
    primer3 = Primer3()
    report = Report()

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


    # Sequence retrieval and markup

    if FUSION:
        # If fusion is passed, run functions required to get the sequences
        # and mark them
        coord_dicts = fusion.split_input(fusion_coords)
        seqs = fusion.fetch_seqs(coord_dicts)

        for ids, dicts in seqs.items():
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
        target_sequence = sequence.fetch_region(
            chrom, startpos - FLANK, endpos + FLANK
        )
        marked_sequence, tagged_string = sequence.markup_sequence(
            FLANK, target_sequence, FUSION, chrom, startpos, endpos
        )

        if startpos == endpos:
            region_id = "{}_{}".format(chrom, startpos)
        else:
            region_id = "{}_{}_{}".format(chrom, startpos, endpos)


    # Creating primers, choosing the best ones and converting to strings
    primer3_results = primer3.run_primer3(region_id, marked_sequence)

    passed_primers = primer3.check_primers(
        region_id, target_sequence, primer3_results, chrom,
        startpos, endpos, FUSION, seqs
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
    output_dir = f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/output/'

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

    for filename in TMP_FILES:
        # clean up temporary files
        print("deleting tmp file: {}".format(filename))
        os.remove(filename)


if __name__ == '__main__':
    main()
