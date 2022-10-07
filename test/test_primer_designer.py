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
    args.total = 1 if args.total == None else args.total
    random_sample = genes_df.sample(n=args.total)

    # list to gather errors in and write to file
    errors = []
    logs = []

    # path to output directory
    test_out = f'{os.path.dirname(os.path.abspath(__file__))}/test_output/'
    print(f'Output Directory: {test_out}')

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
            "-v /home/jason/github/reference:/reference "
            f"-v /home/jason/github/primer_designer/test/test_output:/home/primer_designer/output "
            "--env-file /home/jason/github/primer_designer/.env "
            "primer_designer:2.0.2 "
            f"python -u bin/primer_designer_region.py --chr {chr} --pos {random_pos} --grch38"
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
            print('STDERR', stderr)

        counter += 1

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

