"""
Script to call primer_design_region.py to generate primer designs.
Outputs a zip file of designs.

Kim Brugger (11 Jul 2016), contact: kim@brugger.dk
Nikita Povarnitsyn (30 Jul 2019)
Jethro Rainford (23 Mar 2021)
"""
import os
from pathlib import Path
import random
import re
import string
import subprocess
import sys


def random_string(len=10):
    """
    Creates a random string, used to name out dir uniquely

    Args:
      len (int): length of string to generate, default 10

    returns:
        random_string (str): random character string
    """
    random_string = ''.join(random.choices(
        string.ascii_lowercase + string.ascii_uppercase + string.digits, k=len
    ))

    return random_string


if __name__ == '__main__':

    if len(sys.argv) == 1:
        print("USAGE: bulk_design.py input-file [working dir]")
        print((
            "input-file one tab seperated output-name, region and reference "
            "(grch37 or grch38) of interest pr line"
        ))
        sys.exit()

    infile = sys.argv[1]
    if infile.find(".txt") == -1:
        print("Infile '{}' does not look like a txt file".format(infile))
        sys.exit()

    if len(sys.argv) == 3:
        os.chdir(sys.argv[2])

    outfile = re.sub(r'.txt', '.zip', infile)

    working_dir = random_string()
    os.mkdir(working_dir)
    os.chdir(working_dir)

    primer_designer_path = \
        f"{Path(__file__).parent.absolute()}/primer_designer_region.py"

    with open("../{}".format(infile), "rU") as plist:
        # loop over input txt file of regions to design primers for and
        # call primer_designer_region for each
        for line in plist:
            line = line.strip("\n")
            if line == '':
                continue

            testname, region, reference = re.split(r"[\t ]+", line)
            reference = reference.lower()
            chrom, npos = region.split(":")

            cmd = (
                f"{primer_designer_path} -c {chrom} -p {npos} -o {testname} "
                f"--{reference}"
            )

            try:
                # call primer designer with formatted command
                subprocess.check_call(cmd, shell=True)
            except subprocess.CalledProcessError as call:
                print((
                    f"System call '{cmd}' failed. Exit code: {call.returncode}"
                    f" Error output: {call.output}"
                ))

    os.chdir("..")
    subprocess.check_call("zip -j {} {}/*pdf".format(outfile, working_dir))

    print("SUCCESS\nOutput file: {}".format(outfile))
