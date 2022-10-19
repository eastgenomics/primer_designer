import re
import os
import subprocess
import errno
from pathlib import Path
from shutil import which

from constants import ALLOWED_MISMATCHES, FLANK


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

    def __init__(self, reference: str):
        # get thermoparams dir from primer3 build location
        try:
            path = Path(which("primer3_core")).parents[0]
            self.thermo_params = f'{path}/primer3_config/'
        except TypeError:
            # will raise TypeError from Path if not valid
            print('primer3 missing from path.')

        self.reference = reference

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
            PRIMER_THERMODYNAMIC_PARAMETERS_PATH={self.thermo_params}
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

    def run_primer3(self, seq_id, seq, tmp_files, primer3_file=""):
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
        tmp_files.append(primer3_file)
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

    def check_primers(
            self,
            region_id,
            target_region,
            primer3_dict,
            chromo,
            startpos,
            endpos,
            tmp_files,
            fusion=False,
            seqs=None) -> dict:
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
        tmp_files.append(primers_file)

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
        tmp_files.append(smalt_out)

        # .sma index for reference generated by smalt sometimes has .fasta
        # extension removed, check for both and raise error if can't find one
        ref_file = os.path.expanduser(self.reference).rstrip(".fasta")

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

        if fusion:
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
