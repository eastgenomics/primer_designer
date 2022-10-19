import subprocess
import re

from constants import FLANK, TARGET_LEAD


class Sequence():
    """
    Sequence related functions

    - markup_sequence(): annotates seq., includes marking repeats and snps
    - markup_repeats(): marks repeats longer than 6 bases
    - markup_SNPs(): marks known SNPs in given dbSNP file
    - fetch_known_SNPs(): retrieves kown SNPs from dbSNP file for given region
    """
    def __init__(self, reference, snp):
        self.reference = reference
        self.snp = snp

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

            dbSNPs = self.fetch_known_SNPs(self.snp, chrom, startpos, endpos)
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
                self.snp, chrom, startpos - flank, endpos + flank
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
