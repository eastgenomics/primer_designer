import re
import subprocess
from datetime import datetime
from pathlib import Path

from reportlab.pdfbase.pdfmetrics import stringWidth

from constants import COLOURS, FLANK


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

    def __init__(
            self,
            fusion: bool,
            reference: str,
            primer_text: str,
            primer_version: str,
            snp: str):
        self.fusion = fusion
        self.reference = reference
        self.primer_text = primer_text
        self.primer_version = primer_version
        self.snp = snp

    def make_primer_mapped_str(
            self,
            target_sequence: str,
            passed_primer_seqs: dict):
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

            while (1):
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

    def check_if_primer_clash(self, mapped_list: list, start: int, end: int):
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

    def align_primers_to_seq(self, seq: str, all_primers: list):
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

    def revDNA(self, string: str):
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
            self,
            c,
            y_offset: int,
            primer3_results,
            passed_primers,
            width,
            args,
            FUSION=False,
            chrom=None,
            startpos=None,
            endpos=None,
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
                COLOURS[int(primer_nr)][0], COLOURS[int(primer_nr)][1],
                COLOURS[int(primer_nr)][2]
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

                if (
                    coord_dict[0]['SIDE'] == "<" and
                    coord_dict[0]['STRAND'] == "-1"
                        ):
                    # the coordinates are printed FLANK bases down the given
                    # position s
                    base1 = int(coord_dict[0]['POS']) + FLANK

                elif (
                    coord_dict[0]['SIDE'] == ">" and
                    coord_dict[0]['STRAND'] == "1"
                        ):
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

                if (
                    coord_dict[1]['SIDE'] == "<" and
                    coord_dict[1]['STRAND'] == "-1"
                        ):
                    # the coordinates are printed FLANK bases down the given
                    # position s
                    base1 = int(coord_dict[1]['POS']) + FLANK

                elif (
                    coord_dict[1]['SIDE'] == ">" and
                    coord_dict[1]['STRAND'] == "1"
                        ):
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
            if self.fusion:
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
                            COLOURS[primer_colour[k]][0],
                            COLOURS[primer_colour[k]][1],
                            COLOURS[primer_colour[k]][2]
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
        Function to write common snp annotation text in pdf

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

        if self.fusion:
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
        cmd_args = [
            '; '.join(cmd_args[i:i + 3]) for i in range(0, len(cmd_args), 3)]

        lines.append('')
        lines.append(
            f'Created at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}'
        )
        lines.append(f'SNP file used: {Path(self.snp).stem}')

        if len(Path(self.reference).stem) > 70:
            # break long string into chunks of 70
            long_string = Path(self.reference).stem
            texts = [
                long_string[i:i+70] for i in range(0, len(long_string), 70)]
            for idx, val in enumerate(texts):
                if idx == 0:
                    # first line
                    lines.append(
                        f'Human reference version: {ref} ({val}')
                elif idx == len(texts) - 1:
                    # last line
                    lines.append(f'{val})')
                else:
                    # middle lines
                    lines.append(val)
        else:
            lines.append(
                'Human reference version: {} ({})'.format(
                    ref, Path(self.reference).stem))

        # add final footer text, split over lines as line breaks refused
        # to work and I don't care enough to debug reportlab
        lines.append('')
        for line in self.primer_text:
            lines.append(line)

        lines.append('')
        lines.append(f'primer-designer version: {self.primer_version}')
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
                passed_primers["PRIMER_" + name + "_SEQUENCE"][
                    'MAPPING_SUMMARY']
            ]))

        lines.append("\n")
        lines.append("Consensus sequence:\n")
        lines. append(target_sequence)
        fh.write("\n".join(lines))
        fh.close()
