import re

from constants import FLANK
from utils import fetch_region


class Fusion():
    """
    Fusion primer designing related functions

    - split_input(): creates dicts for given coords, used to add attributes
    - fetch_seqs(): calls Sequence.fetch_region() to get sequence from ref
    - flip_fusion_seq(): returns reverse complement where seq is on same
        strand either side of breakpoint
    """
    def __init__(self, reference):
        self.reference = reference

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
                coord_dict[index]['SEQ'] = fetch_region(
                    dict['CHR'],
                    dict['POS'],
                    dict['POS'] + FLANK,
                    self.reference
                )
            elif dict['SIDE'] == ">":
                coord_dict[index]['SEQ'] = fetch_region(
                    dict['CHR'],
                    dict['POS'] - FLANK,
                    dict['POS'],
                    self.reference
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

            elif (
                seqs_dict[i]['SIDE'] == ">" and
                seqs_dict[i]['STRAND'] == "-1"
                    ):
                # on rev strand, get reverse complement of nroaml and marked
                # sequence, reverse tagged seq and add DARK_SIDE
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])

                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])

                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = "<"

            elif (
                seqs_dict[i]['SIDE'] == "<" and
                seqs_dict[i]['STRAND'] == "-1"
                    ):
                # on rev strand, get reverse complement of nroaml and marked
                # sequence, reverse tagged seq and add DARK_SIDE
                seqs_dict[i]['SEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['SEQ'][::-1])

                seqs_dict[i]['MSEQ'] = "".join(complement.get(
                    base, base) for base in seqs_dict[i]['MSEQ'][::-1])

                seqs_dict[i]['TSEQ'] = seqs_dict[i]['TSEQ'][::-1]
                seqs_dict[i]['DARK_SIDE'] = ">"

        if (
            seqs_dict[0]['DARK_SIDE'] == ">" and
            seqs_dict[1]['DARK_SIDE'] == "<"
                ):
            # one set seqs before and after breakpoint => valid (?)
            target_sequence = ''.join(
                [seqs_dict[0]['SEQ'], seqs_dict[1]['SEQ']])
            marked_sequence = ''.join(
                [seqs_dict[0]['MSEQ'], seqs_dict[1]['MSEQ']])
            tagged_string = ''.join(
                [seqs_dict[0]['TSEQ'], seqs_dict[1]['TSEQ']])

        elif (
            seqs_dict[0]['DARK_SIDE'] == "<" and
            seqs_dict[1]['DARK_SIDE'] == ">"
                ):
            # one set seqs before and after breakpoint => valid (?)
            target_sequence = ''.join(
                [seqs_dict[1]['SEQ'], seqs_dict[0]['SEQ']])
            marked_sequence = ''.join(
                [seqs_dict[1]['MSEQ'], seqs_dict[0]['MSEQ']])
            tagged_string = ''.join(
                [seqs_dict[1]['TSEQ'], seqs_dict[0]['TSEQ']])

        else:
            # I think it gets here if both of the coordinates are either
            # before or after the breakpoint?
            raise ValueError((
                "An error occured in designing fusion primers, this design "
                "does not seem possible as both sequences are on the same "
                "side of the breakpoint. Please review the parameters passed."
            ))

        return target_sequence, tagged_string, marked_sequence
