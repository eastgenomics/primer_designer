import re
import subprocess


def fetch_region(chrom: str, start: int, end: int, reference: str):
    """
    Gets the nucleotide sequence from the appropriate reference
    file with the given coordinates

    Args:
        - chrom (str): target chromosome
        - start (int): target start coord
        - end (int): target end coord
        - reference (str): path for ref file

    Returns:
        - sequence (str): target sequence from reference file

    """
    # call samtools to get sequence from fasta
    cmd = f"samtools faidx {reference}  {chrom}:{start}-{end}"

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
