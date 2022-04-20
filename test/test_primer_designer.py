"""
Script to test primer designer in a quick and hacky way by generating lots
of reports, if any fail traceback will be wrriten out to file.
"""
import argparse
import glob
import os
import pandas as pd
import random
import subprocess
import sys
import traceback


def load_bed():
    """
    Load in bed file to generate positions from

    """
    bed_file = f'{os.path.dirname(os.path.abspath(__file__))}/test_genes.bed'

    with open(bed_file) as bed:
        genes_df = pd.read_csv(
            bed, sep='\t', names=['chr', 'start', 'end', 'tx']
        )

    return genes_df


def run_test(args, genes_df):
    """
    Gets rows randomly from bed file df, generates random position between
    the start and end then runs primer designer

    Args:
        - genes_df (df): df of bed file
    """
    # get random sample of rows
    random_sample = genes_df.sample(n=args.total)

    # list to gather errors in and write to file
    errors = []
    logs = []

    # path to output directory
    test_out = f'{os.path.dirname(os.path.abspath(__file__))}/test_output/'

    counter = 1
    for _, row in random_sample.iterrows():
        # generate random position between coordinates of row
        chr = row['chr']
        random_pos = random.randint(row['start'], row['end'])

        print(
            f'Generating report {counter}/{args.total} for {chr}:{random_pos}'
        )
        logs.append(
            f'Generating report {counter}/{args.total} for {chr}:{random_pos}'
        )

        cmd = (
            "docker run "
            "-v /home/$USER/reference_files:/reference_files "
            f"-v {test_out}:/home/primer_designer/output "
            "--env REF_37=/reference_files/grch37/human_g1k_v37.fasta "
            "--env DBSNP_37=/reference_files/grch37/dbsnp_grch37_common_all_20180423.vcf.gz "
            "primer_designer "
            f"python3 bin/primer_designer_region.py -c {chr} -p {random_pos} --grch37"
        )

        # call primer designer docker to generate report
        try:
            subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT
            )

            if args.delete:
                # delete last generated pdf
                os.remove(
                    max(glob.glob(rf'{test_out}*.pdf'), key=os.path.getctime)
                )
        except subprocess.CalledProcessError as exc:
            print('Oh no, an error, adding to file')
            errors.append(f"Error designing primers for {chr}:{random_pos}")
            stderr = exc.output.decode()
            errors.append(stderr)

        counter += 1

    with open('test.log', 'w') as log_file:
        # write a log of everything generated
        for log in logs:
            log_file.write(str(log))
            log_file.write('\n')

    if len(errors) > 0:
        # found some errors, write to file
        with open('test_errors.txt', 'w') as err_file:
            for error in errors:
                for line in error:
                    err_file.write(str(line))
                err_file.write('\n')


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-t', '--total', type=int,
        help='Total number of reports to generate'
    )
    parser.add_argument(
        '--delete', action='store_true',
        help='delete generate PDF after it is created'
    )

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    genes_df = load_bed()
    run_test(args, genes_df)

