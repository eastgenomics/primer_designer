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
from config import SAMTOOLS, TABIX, SMALT, PRIMER3, FONT, THERMO_PARAMS,\
    REF_37, REF_38, DBSNP_37, DBSNP_38

# for config
THERMO_PARAMS = "/mnt/storage/apps/software/primer3_core/2.3.7/src/primer3_config/"

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

        Args:
            - flank (int):
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

        return sequence_list, tagged_string


    def markup_repeats(self, tags, sequence):
        """
        Marks the repeats in the sequence longer than 6 bases.

        Args:
            - tags (list): empty list to add tags to (?)
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
                if re.search('<', sequence_list[mask_pos + i]) or re.search(
                    r'\[', sequence_list[mask_pos + i]
                ):
                    continue

                sequence_list[mask_pos + i] = ' <' + sequence_list[mask_pos + i] + '> '
                tags[mask_pos + i] = 'X'
            else:
                pass

        sequence_list = "".join(sequence_list)
        sequence_list = re.sub(' ', '', sequence_list)
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
        cmd = f"{TABIX} {tabix_file}  {chrom}:{start}-{end}"

        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)

        output = p.communicate()
        snps = []

        for line in output[0].split("\n"):
            snps.append(line.split("\t"))

        return snps


class Primer3():
    """Functions for calling & handling output from primer3"""

    def template(self, seq_id, seq, flank, len):
        """
        Returns str of template params for generating pdf, in separate
        function as it is long and ugly.

        Args:
            - seq_id (str):
            - seq (str):
            - flank (int):
            - len (int):
        Returns:
            - template (str): formatted str of template parameters
        """
        template = f'''
            SEQUENCE_ID={seq_id}
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
            =
        '''

        return template


    def run_primer3(self, seq_id, seq, primer3_file=""):
        """
        Creates a primer3 file. Runs PRIMER and captures its output and returns it as a dictionary where the name of the
        primer is the key and the sequence is the value
        """
        verbose_print("run_primer3", 2)

        if primer3_file == "":
            primer3_file = seq_id + ".primer3"

        primer3_file = re.sub("[<>:]", "_", primer3_file)

        TMP_FILES.append(primer3_file)

        write_primer3_file(seq_id, seq, primer3_file)

        cmd = "{} < {}".format(PRIMER3, primer3_file)

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        output = process.communicate()
        output_dict = dict()

        for line in output[0].split("\n"):
            if line == '=':
                break

            key, value = line.split("=")

            output_dict[key] = value

        return output_dict


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
        len = downstream_start - upstream_end + 1
        template = self.template(seq_id, seq, flank, len)

        with open(primer3_file, 'w+') as outfile:
            outfile.write(template)

        outfile.close()


    def check_primers(self, region_id, target_region, primer3_dict,
                      chromo, startpos, endpos, FUSION=False, seqs=None):
        """
        Creates a fasta file with sequences of the primers, which are
        requrired by SMALT to map to the reference genome. The smalt result
        is saved in a nested dictionary containing mapping summary of each
        primer.

        Args:
            - region_id (int):
            - target_region (str):
            - primer3_dict (dict):
            - chrom (str):
            - startpos (int):
            - endpos (int):
            - FUSION (bool):
            - seqs (str)
        Returns:
            - seq_dict (dict):
        """
        verbose_print("check_primers", 2)

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

                primerfasta.write(">{}\n{}\n".format(pl_id, primer3_dict[pl_id]))
                primerfasta.write(">{}\n{}\n".format(pr_id, primer3_dict[pr_id]))

            primerfasta.close()

        smalt_out = region_id + ".smalt"
        TMP_FILES.append(smalt_out)

        ref = re.sub(r'\..*', '', REFERENCE)

        cmd = f"{SMALT} map  -d -1 -m 15 {ref} {primers_file} > {smalt_out}"

        verbose_print(cmd, 4)
        subprocess.call(cmd, shell=True)

        seq_dict = {}

        with open(smalt_results, "r") as f:
            for line in f.readlines():

                if line.startswith("@") or line.startswith("FULLSEQ"):
                    continue

                line = line.split('\t')

                name = line[0]
                chrom = line[2]
                pos = line[3]
                mismatch = len(line[9]) - int(line[12].split(':')[2])

                if mismatch > 5:
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

            for ids, dicts in seqs.items():
                if dicts["DARK_SIDE"] == ">":
                    chromo = dicts['CHR']
                    startpos = int(dicts['POS']) - FLANK
                    endpos = int(dicts['POS'])
                    left_dict = self.primer_report(chromo, startpos, endpos, left_dict)

                elif dicts["DARK_SIDE"] == "<":
                    chromo = dicts['CHR']
                    startpos = int(dicts['POS'])
                    endpos = int(dicts['POS']) + FLANK
                    right_dict = self.primer_report(chromo, startpos, endpos, right_dict)

            # merge two dictionaries back into one
            seq_dict = self.merge_dicts(left_dict, right_dict)

        else:
            # declare chrom ,start ,end
            seq_dict = self.primer_report(chromo, startpos, endpos, seq_dict)

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
        result = {}
        for dictionary in dicts:
            result.update(dictionary)
        return merged_dicts


    def primer_report(self, chromo, startpos, endpos, seq_dict):
        """
        Creating a nested dictionary with summary report for each primer

        Args:
            - chrom (str):
            - startpos (int):
            - endpos (int):
            - seq_dict (dict):
        Returns:
            - seq_dict (dict):
        """
        for k, v in seq_dict.iteritems():
            uniq_chr = []
            mis_chr = []

            for pos, chrom in enumerate(seq_dict[k]['CHR']):

                if seq_dict[k]['CHR'][pos] == chromo and int(
                        seq_dict[k]['POS'][pos]) in range(
                        startpos - FLANK,
                        endpos + 1 + FLANK) and seq_dict[k]['MISMATCH'][pos] == '0':
                    uniq_chr.append(seq_dict[k]['CHR'][pos])

                else:
                    mis_chr.append(seq_dict[k]['CHR'][pos])

            if len(uniq_chr) == 1 and len(mis_chr) == 0:

                seq_dict[k]['MAPPING_SUMMARY'] = 'unique mapping'

            elif len(uniq_chr) + len(mis_chr) in range(2, 6):
                seq_dict[k]['MAPPING_SUMMARY'] = '{} mapping(s) on chr {} with {} mismatches'.format(
                    len(uniq_chr + mis_chr), ', '.join(uniq_chr + mis_chr), ', '.join(mydict[k]['MISMATCH']))

            elif len(uniq_chr) + len(mis_chr) > 5:
                seq_dict[k]['MAPPING_SUMMARY'] = '{} mapping(s)'.format(
                    len(uniq_chr + mis_chr))

        return seq_dict


    def extract_passed_primer_seqs(self, primer3_results, passed_primers):
        """
        Extracting the name and sequnece of the primers generated by primer3

        Args:
            - primer3_results ():
            - passed_primers ():
        Returns:
            - primer_seqs (list):
        """
        verbose_print("extract_passed_primer_seqs", 2)

        primer_seqs = []
        left_primer_seqs = []
        right_primer_seqs = []
        for primer in passed_primers:

            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            primer_seqs.append([name, primer3_results[primer]])

        return primer_seqs


    def make_primer_mapped_strings(self, target_sequence, passed_primer_seqs):
        """
        Creates a mapping string out of the target sequence and primers.
        Used to be printed on the PDF report.

        Args:
            - target_sequence ():
            - passed_primer_seqs ():
        Returns:
            - mapped_strings (list):
            - mapped_colours (list):
        """
        mappings = align_primers_to_seq(target_sequence, passed_primer_seqs)

        mapped_strings = []
        mapped_strings.append([" "] * len(target_sequence))

        mapped_colours = []
        mapped_colours.append([-1] * len(target_sequence))

        for mapping in mappings:
            name, primer, pos, strand = mapping

            primer_nr = int(re.sub(r'.*_(\d)', r'\1', name))
            mapping_index = 0

            arrows = (len(primer) - len(name) - 2) / 2
            arrow_type = ">"

            if strand == 1:
                # if the primer is on the reverse strand
                primer = revDNA(primer)
                arrow_type = "<"

            tag = arrow_type * arrows + " " + name + " " + arrow_type * arrows

            if len(tag) < len(primer):
                tag += arrow_type

            primer = list(tag)
            scan_index = 0

            while(1):
                if check_if_primer_clash(
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


    def align_primers_to_seq(self, seq, all_primers):
        """
        Checking if the primers and their revers sequences match any
        regions of the sequence of interest adding everything to a list.
        Used in make_primers_mapped_strings().

        Args:
            - seq ():
            - all_primers():
        Returns:
            - mappings ():
        """
        verbose_print("align_primers_to_seq", 2)
        primers = []
        rev_primers = []

        for primer_set in all_primers:
            name, primer = primer_set
            primers.append([name, primer])
            rev_primers.append([name, revDNA(primer)])

        mappings = []

        for i in range(0, len(seq)):

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
        verbose_print("check_if_primers_clash", 4)

        for i in range(start, end):
            if mapped_list[i] != " ":
                return True

        return False


    def revDNA(self, string):
        """
        Returns reverse compliment of given nucleotide sequence.

        Args:
            - string (str): nucleotide sequence
        Returns:
            - rev_str (str): nucleotide sequence but magically reversed
        """
        verbose_print("revDNA", 4)
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


class Report():
    """Functions to generate PDF report"""

    def pretty_pdf_primer_data(
            self, c, y_offset, primer3_results, passed_primers, width,
            FUSION=False, chrom=None, startpos=None, endpos=None,
            coord_dict=None):
        """
        Adds the header primer information to the output PDF file

        Args:
            -
        Returns:
            -
        """
        verbose_print("pretty_pdf_primer_data", 2)
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8

        if FUSION:
            c.drawString(
                40,
                y_offset,
                "Primer design report for a fusion between chr: {} position: {} and chr: {} position: {} ".format(
                    coord_dict[0]['CHR'],
                    coord_dict[0]['POS'],
                    coord_dict[1]['CHR'],
                    coord_dict[1]['POS']))

        else:
            if startpos == endpos:
                c.drawString(
                    40,
                    y_offset,
                    "Primer design report for chr {} position: {}".format(
                        chrom,
                        startpos))
            else:
                c.drawString(
                    40, y_offset,
                    "Primer design report for chr: {} range: {}-{}".format(
                        chrom, startpos, endpos
                    )
                )

        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 16

        c.drawString(
            40, y_offset,
            "ID         %GC    TM     Primer sequence                       Mapping(s)")
        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8

        primer_seqs = []

        for primer in sorted(passed_primers):
            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            if name == "RIGHT_0":
                y_offset -= 8

            picked_primer = ' '

            c.drawString(40,
                        y_offset,
                        "{:10} {:.2f}  {:.2f}  {:25}             {}            ".format("",
                                                                                        float(primer3_results["PRIMER_" + name + "_GC_PERCENT"]),
                                                                                        float(primer3_results["PRIMER_" + name + "_TM"]),
                                                                                        primer3_results["PRIMER_" +name + "_SEQUENCE"],
                                                                                        passed_primers["PRIMER_" + name +"_SEQUENCE"]['MAPPING_SUMMARY'],
                                                                                        ))

            primer_nr = re.sub(r'.*_(\d)', r'\1', name)

            c.setFillColorRGB(colours[int(primer_nr)][0],
                            colours[int(primer_nr)][1],
                            colours[int(primer_nr)][2])

            c.drawString(40, y_offset, name)
            c.setFillColorRGB(0, 0, 0)
            y_offset -= 8

        y_offset -= 8
        c.line(40, y_offset, width - 40, y_offset + 2)
        y_offset -= 8
        y_offset -= 8

        return y_offset


    def pretty_pdf_fusion_mappings(self, top_offset, c, coord_dict, passed_primer_seqs, FUSION=False):
        """Function to interogate the nested dictionary. Firstly the first sequence with the
        position of interest at the end is printed. It is followed by the other one. DARK_SIDE
        varible allows to track the change of the side of the flipped sequence. If sequence of interest
        was after given position and then reverse-complemented the signs would be < > for side and dark_side
        respectively.

        """

        verbose_print("pretty_pdf_fusion_mappings", 2)
        # spaces required to match the positions of two breakpoints when pdf is
        # created
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
                primer_strings, primer_colours = make_primer_mapped_strings(
                    target_sequence, passed_primer_seqs)
                top_offset = pretty_pdf_mappings(
                    top_offset,
                    target_sequence,
                    tagged_string,
                    primer_strings,
                    primer_colours,
                    base1,
                    c,
                    side,
                    darkside)

                target_sequence = spaces + coord_dict[1]['SEQ']
                tagged_string = spaces + coord_dict[1]['TSEQ']
                base1 = int(coord_dict[1]['POS'])
                side = coord_dict[1]['DARK_SIDE']
                darkside = coord_dict[1]['DARK_SIDE']
                primer_strings, primer_colours = make_primer_mapped_strings(
                    target_sequence, passed_primer_seqs)
                top_offset = pretty_pdf_mappings(
                    top_offset,
                    target_sequence,
                    tagged_string,
                    primer_strings,
                    primer_colours,
                    base1,
                    c,
                    side,
                    darkside)

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
                primer_strings, primer_colours = make_primer_mapped_strings(
                    target_sequence, passed_primer_seqs)
                top_offset = pretty_pdf_mappings(
                    top_offset,
                    target_sequence,
                    tagged_string,
                    primer_strings,
                    primer_colours,
                    base1,
                    c,
                    side,
                    darkside)

                target_sequence = spaces + coord_dict[0]['SEQ']
                tagged_string = spaces + coord_dict[0]['TSEQ']
                base1 = int(coord_dict[0]['POS'])
                side = coord_dict[0]['DARK_SIDE']
                darkside = coord_dict[0]['DARK_SIDE']
                primer_strings, primer_colours = make_primer_mapped_strings(
                    target_sequence, passed_primer_seqs)
                top_offset = pretty_pdf_mappings(
                    top_offset,
                    target_sequence,
                    tagged_string,
                    primer_strings,
                    primer_colours,
                    base1,
                    c,
                    side,
                    darkside)

        return (top_offset)


    def pretty_pdf_mappings(
            self,
            top_offset,
            target_sequence,
            tagged_string,
            primer_strings,
            primer_colours,
            base1,
            c,
            side=None,
            darkside=None):
        """
        Outputs the mapping results into the .pdf file
        """

        verbose_print("pretty_pdf_mappings", 2)
        spaces = ' ' * ((FLANK + 1) % 80)

        # The flip of one of the sequences defines how the numeration of
        # coordinates go, either in descending or ascending orders

        for i in range(0, len(target_sequence), 80):

            if FUSION:

                if side == '>' and darkside == '>':
                    p_line = "{}  {}".format(base1 + i, target_sequence[i: i + 80])
                    m_line = "           " + tagged_string[i: i + 80]

                elif side == '<' and darkside == '<':
                    p_line = "{}  {}".format(
                        base1 + i - len(spaces), target_sequence[i: i + 80])
                    m_line = "           " + tagged_string[i: i + 80]

                elif side == '>' and darkside == '<':
                    p_line = "{}  {}".format(
                        base1 - i + len(spaces), target_sequence[i: i + 80])
                    m_line = "           " + tagged_string[i: i + 80]

                elif side == '<' and darkside == '>':
                    p_line = "{}  {}".format(base1 - i, target_sequence[i: i + 80])
                    m_line = "           " + tagged_string[i: i + 80]

            else:
                p_line = "{}  {}".format(base1 + i, target_sequence[i: i + 80])
                m_line = "           " + tagged_string[i: i + 80]

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

                c.drawString(x_offset, top_offset, p_line[k])
                x_offset += stringWidth(" ", 'mono', 8)
                c.setFillColorRGB(0, 0, 0)

            top_offset -= 8

            m_line = re.sub(r'X', r' ', m_line)

            if re.search(r'\*', m_line) or re.search(r'\-', m_line):
                # setting colour for the TARGET LEAD (50 base pairs around the
                # position of interest)
                c.setFillColorRGB(0, 190, 0)
                c.drawString(40, top_offset, m_line)
                c.setFillColorRGB(0, 0, 0)

                top_offset -= 8

            for j in range(0, len(primer_strings)):
                primer_string = primer_strings[j]
                primer_colour = primer_colours[j]

                line = primer_string[i: i + 80]
                if (re.match(r'^ *$', line)):
                    continue

                x_offset = 40 + stringWidth(" ", 'mono', 8) * 11

                for k in range(i, i + 80):
                    if k > len(target_sequence) - 1:
                        break

                    if primer_colour[k] >= 0:

                        c.setFillColorRGB(colours[primer_colour[k]][0],
                                        colours[primer_colour[k]][1],
                                        colours[primer_colour[k]][2])

                    c.drawString(x_offset, top_offset, primer_string[k])
                    x_offset += stringWidth(" ", 'mono', 8)
                    c.setFillColorRGB(0, 0, 0)

                top_offset -= 8

            top_offset -= 8

        return (top_offset)


    def pretty_pdf_method(self, top_offset, args, c, seqs=None):
        verbose_print("pretty_pdf_method", 3)

        lines = method_blurb(args, seqs)

        top_offset = 80

        for line in lines:
            c.drawString(40, top_offset, line)
            top_offset -= 8


    def method_blurb(args, seqs=None):

        lines = []

        if FUSION:

            lines.append(
                'Primer design report for a fusion between chromosome ' +
                seqs[0]['CHR'] +
                ' at the position ' +
                seqs[0]['POS'])
            lines.append(
                'and chromosome ' +
                seqs[1]['CHR'] +
                ' at the position ' +
                seqs[0]['POS'] +
                '.')
            for regionid, region_dict in seqs.items():

                if region_dict['STRAND'] == "-1":
                    reverse_region = region_dict['CHR']
                    lines.append(
                        'The sequence on chromosome ' +
                        reverse_region +
                        ' was reverse-complemented to produce this fusion sequence')

                else:
                    continue

        if args.grch37 or args.hg19:
            lines.append(
                'primer-designer version: ' +
                VERSION +
                ' using dbSNP 144 for SNP checking, and human reference GRCh37.')
        else:
            lines.append(
                'primer-designer version: ' +
                VERSION +
                ' using dbSNP 150 for SNP checking, and human reference GRCh38.')

        lines.append(
            'Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
        lines.append(
            'allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

        return (lines)


    def pretty_primer_data(
            outfile,
            primer3_results,
            passed_primers,
            chrom,
            startpos,
            endpos,
            target_sequence=None,
            FUSION=False,
            coord_dict=None):
        """
        function to output the report in txt format. No primers would be
        displayed on the target sequence
        """
        verbose_print("pretty_primer_data", 2)

        fh = open(outfile, 'w')

        lines = []
        if FUSION:
            lines.append(
                "Primer design report for a fusion between chr: {} position: {} and chr: {} position: {} ".format(
                    coord_dict[0]['CHR'],
                    coord_dict[0]['POS'],
                    coord_dict[1]['CHR'],
                    coord_dict[1]['POS']))

        elif startpos == endpos:
            lines.append(
                "Primer design report for chr: {} position: {}".format(
                    chr, startpos))

        else:
            lines.append(
                "Primer design report for chr: {} range: {}-{}".format(chr, startpos, endpos))

        lines.append("ID\t%GC\tTM\tPrimer sequence\tBest primer\tMapping(s)")

        for primer in sorted(passed_primers):
            if primer == 'FULLSEQ':
                continue

            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            lines.append("\t".join([name, primer3_results["PRIMER_" + name + "_GC_PERCENT"],
                                    primer3_results["PRIMER_" + name + "_TM"],
                                    primer3_results["PRIMER_" +
                                                    name + "_SEQUENCE"],
                                    passed_primers["PRIMER_" + name + "_SEQUENCE"]['MAPPING_SUMMARY']]))

        lines.append("\n")
        lines.append("Consensus sequence:\n")
        lines. append(target_sequence)
        fh.write("\n".join(lines))
        fh.close()


    def pretty_print_mappings(target_sequence, tagged_string, primer_strings, base1):
        """used to print out the information of the primers in the terminal"""

        lines = []

        for i in range(0, len(target_sequence), 80):

            lines.append("{}  {}".format(base1 + i, target_sequence[i: i + 80]))
            lines.append("           " + tagged_string[i: i + 80])

            for primer_string in primer_strings:

                line = "           " + primer_string[i: i + 80]

                if re.match(r'^ *$', line):
                    continue

                lines.append(line)

            lines.append("")

        return (lines)


    def pretty_print_primer_data(primer3_results, passed_primers):
        """
        function which is never used.
        """

        lines = []

        verbose_print("extract_passed_primer_seqs", 3)

        lines.append("\n")
        lines.append("\n")
        lines.append("\n")
        lines.append("_-=-" * 15 + "_")

        if startpos == endpos:
            lines.append(
                " Primer design report for chr: {} position: {}".format(
                    chr, startpos))
        else:
            lines.append(
                " Primer design report for chr: {} range: {}-{}".format(chr, startpos, endpos))

        lines.append("_-=-" * 15 + "_")
        lines.append("\n")

        lines.append(
            "ID         %GC    TM     Primer sequence           Mapping(s)    ")
        lines.append("_-=-" * 15 + "_")

        primer_seqs = []
        for primer in sorted(passed_primers):
            if primer == 'FULLSEQ':
                continue

            name = primer
            name = re.sub(r'PRIMER_', '', name)
            name = re.sub(r'_SEQUENCE', '', name)

            lines.append("{} {}  {}  {} {}".format(name,
                                                float(primer3_results["PRIMER_" + name + "_GC_PERCENT"]),
                                                float(primer3_results["PRIMER_" + name + "_TM"]),
                                                primer3_results["PRIMER_" +name + "_SEQUENCE"],
                                                passed_primers["PRIMER_" + name + "_SEQUENCE"]['MAPPING_SUMMARY']))

        lines.append("")
        lines.append("-=" * 46)
        lines.append("")


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
    primer3 = Primer3()
    report = Report()

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

        target_sequence, tagged_string, marked_sequence = \
            fusion.flip_fusion_seq(seqs)

        region_id = fusion
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


    # Creating primers and choosing the best ones; and converting to strings
    primer3_results = run_primer3(region_id, marked_sequence)

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
            primer3.make_primer_mapped_strings(
                target_sequence, passed_primer_seqs
            )

        lines = report.pretty_print_mappings(
            target_sequence, tagged_string, mapped_primer_strings,
            startpos - FLANK
        )

    # generate output PDF report
    filename = re.sub("[<>:]", "_", region_id)

    if args.output:
        filename = args.output 

    if args.text_output:
        report.pretty_primer_data(
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
        top_offset = report.pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width,
            FUSION, chrom, startpos, endpos, seqs
        )

        c.setFont('mono', 8)
        report.pretty_pdf_fusion_mappings(
            top_offset, c, seqs, passed_primer_seqs, FUSION)

        report.pretty_pdf_method(top_offset, args, c, seqs)

    else:
        # normal primer report formatting
        c.setFont('mono', 6)
        top_offset = report.pretty_pdf_primer_data(
            c, height - 30, primer3_results, passed_primers, width, FUSION,
            chrom, startpos, endpos
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
        print("deleting tmp file: {}".format(filename))
        os.remove(filename)


if __name__ == '__main__':
    main()
