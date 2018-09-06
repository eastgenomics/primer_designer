#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (15 Sep 2015), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

import shlex
import subprocess
import re

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
width, height = A4 #keep for later

colours = [[255,   0,   0], # red
           [  0, 255,   0], # green
           [  0,   0, 255], # blue
           [255,   0, 255], # Pink
           [0,   255, 255],
           [255, 255,   0]] # Yellow, crap!

# External tools we use 
SAMTOOLS   = '/software/bin/samtools ';
REFERENCE  = '/refs/human_1kg/human_g1k_v37.fasta '
TABIX      = '/software/bin/tabix '
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '

# default parameters, can be overridden by the user (perhaps)
FLANK              = 500
NR_PRIMERS         = 4
ALLOWED_MISMATCHES = 4
MAX_MAPPINGS       = 5

VERBOSE    =  3
VERSION    =  '1.0-rc2'

def verbose_print( msg, level ):
    if ( level <= VERBOSE ):
        print msg

def fetch_region( chr, start, end ):

    verbose_print( "fetch_region", 3)

    cmd = "%s faidx %s  %s:%d-%d " % ( SAMTOOLS, REFERENCE, chr, start, end )
    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    sequence = ""
    for line in ( output[0].split("\n")):
        if ( re.match('>', line)):
            continue 

        sequence += line
    return sequence 


def fetch_known_SNPs( tabix_file, chr, start, end):
    verbose_print( "fetch_known_SNPs", 3)

    cmd = "%s %s  %s:%d-%d " % ( TABIX, tabix_file, chr, start, end )

    print cmd

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

#    pp.pprint( output );

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 



def write_primer3_file(seq_id, seq, primer3_file = ""):
    verbose_print( "write_primer3_file", 4)

    template = '''SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_TARGET={flank},1
PRIMER_FIRST_BASE_INDEX=1
PRIMER_TASK=generic
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_PRODUCT_SIZE_RANGE=150-500
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
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/software/bin/primer3_config/
='''

    details = {
        "seq_id" : seq_id,
        "seq"    : seq,
        "flank"  : FLANK
    }



    with open( primer3_file, 'w+') as outfile:
        outfile.write(template.format(**details))

    outfile.close()


# This captures the $_ output of the primer3 engine and returns it as a dictionary.
def run_primer3( seq_id, seq, primer3_file = ""):
    verbose_print( "run_primer3", 3)

    if ( primer3_file == "" ):
        primer3_file = seq_id + ".primer3"

    primer3_file = re.sub(":", "_", primer3_file)


    write_primer3_file(seq_id, seq, primer3_file)


    cmd = PRIMER3 + " < " + primer3_file


    print cmd
    args = shlex.split( cmd )
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               shell=True)
  
    output = process.communicate()
    output_dict = dict()
  
    for line in output[0].split("\n"):
    #print line

        if (line == '='):
            break

        key, value = line.split("=")

        output_dict[ key ] = value


  #pp.pprint( output_dict )

    return output_dict



def check_primers( region_id, target_region, primer3_dict):
    verbose_print( "check_primers", 3)


    primers_file = region_id + "_p3seq.fasta"

    with open( primers_file, 'w+') as primerfasta:

        primerfasta.write(">FULLSEQ\n" + target_region + "\n")

        for i in range(0, 5):
            pl_id = "PRIMER_LEFT_%d_SEQUENCE" % i
            pr_id = "PRIMER_RIGHT_%d_SEQUENCE" % i

            if ( pr_id not in primer3_dict or pl_id not in primer3_dict):
                continue

            primerfasta.write(">%s\n%s\n" % (pl_id, primer3_dict[ pl_id ]))
            primerfasta.write(">%s\n%s\n" % (pr_id, primer3_dict[ pr_id ]))

        primerfasta.close()

    smalt_results = region_id + ".smalt"
    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"
    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results 

    print cmd
    subprocess.call(cmd, shell=True)


    id_word = "FULLSEQ"
    match = {}
    smalt_report = []
    query_region = []

    res = dict()

    with open(smalt_results, 'rU') as smalt_output:

        for line in smalt_output:
    
            if (line.startswith("@")):
                continue
       
            line = line.rstrip("\n")
            fields = line.split("\t")
#            pp.pprint( fields )
            match[ 'name'           ] = fields[  0  ] #mapping_score
            match[ 'chromosome'     ] = fields[  2  ] #chromosome
            match[ 'pos'            ] = fields[  3  ] #mapping position
            match[ 'length'         ] = len(fields[  9  ]) #mapping_score
            match[ 'length_matched' ] = int(re.sub("AS:i:", '', fields[  12 ])) #mapping_length


            match_id     = fields[ 0 ] 
            match_chr    = fields[ 2 ]
            match_pos    = fields[ 3 ]
            match_length = len(fields[ 9 ])
            match_mathes = int(re.sub("AS:i:", '', fields[  12 ]))

            if ( match_id not in res ):
                res[ match_id ] = dict()
                res[ match_id ][ 'CHR' ] = []
                res[ match_id ][ 'POS' ] = []



            if (match['length'] <= match['length_matched'] + ALLOWED_MISMATCHES):
                res[ match_id ][ 'CHR' ].append( match_chr )
                res[ match_id ][ 'POS' ].append( match_pos )



#            smalt_report.append()
    


    for primer  in  res.keys() :
        if (primer == 'FULLSEQ'):
            continue 

        res[ primer ]['MAPPING_SUMMARY'] = 'unique mapping'


        nr_of_chromosomes = len(set(res[ primer ][ 'CHR' ]))
        nr_of_mappings    = len( res[ primer ][ 'POS' ])

        if (nr_of_mappings > 1 and nr_of_mappings <= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosomes: " + ",".join(res[ primer ][ 'CHR' ])


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes



#    pp.pprint( smalt_report)
#    pp.pprint( res )


    return res

def revDNA( string ):
    verbose_print( "revDNA", 5)
    rev_bases = { 'A': 'T', 
                  'a': 'T', 
                  'C': 'G',
                  'c': 'G',
                  'G': 'C',
                  'g': 'C',
                  'T': 'A',
                  'T': 'A',
                  '-': '-'}

#    print "STRING FOR REVDNA:: " + string
    rev_str = len(string)*[None]
    for i in range(0, len(string)):
        rev_str[ len(string) - i - 1] =  rev_bases[ string[ i ]]

    return "".join( rev_str )



def align_primers_to_seq( seq, all_primers):
    verbose_print( "align_primers_to_seq", 4)

    primers = []
    rev_primers = []
    for primer_set in all_primers:
        ( name, primer) = primer_set
#        print primer
        primers.append(  [name, primer] )
        rev_primers.append( [name, revDNA( primer )])

#    pp.pprint( primers )

    mappings = []

    for i in range(0, len(seq)):
        for primer_set in primers:
            (name, primer) = primer_set
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
#                print primer + " Matches at pos " + str( i )

#                mappings.append( [primer, i, 0] )
                mappings.append( [name, primer, i, 0] )

        for primer_set in rev_primers:
            (name, primer) = primer_set
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
#                print primer + " Matches at pos " + str( i ) + " minus! "

#                mappings.append( [primer, i, 1] )
                mappings.append( [name, primer, i, 1] )


#    pp.pprint( mappings )

    return mappings



def get_and_parse_arguments():
    verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chr')
    parser.add_argument('-p', '--pos')
    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank')
    parser.add_argument('-t', '--text_output')

    args = parser.parse_args()

    return args



def markup_sequence( chr, pos, flank, sequence):
    verbose_print( "tag_sequence", 3)

    sequence = list(sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    tags[ FLANK ] = '*'

    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ flank ] = ' [' + target_sequence [ flank ] + '] '

    dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/annots-rsIDs-dbSNPv144.20150605.tab.gz', chr, pos - flank, pos + flank )

    masked_positions = []

    for dbSNP in (dbSNPs):
        if ( len( dbSNP ) < 6):
            continue
    
        (snp_chr, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP


        snp_pos = int( snp_pos )
        if ( common == '1'):
            #        pp.pprint( dbSNP )
            mask_pos = snp_pos - (pos - flank) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
            if ( mask_pos == flank ):
                continue
        
        # In the odd case we are looking at common deletion, mask the whole region. Normally this will just be one base
            for i in range(0, len(snp_ref)):
                if (len(sequence)<=  mask_pos + i):
                    break

                # If already masked skip masking it again.
                if ( re.search('<', sequence[ mask_pos + i  ]) or 
                     re.search('\[',sequence[ mask_pos + i  ])):
                    continue

                sequence[ mask_pos + i  ] = ' <' + sequence[ mask_pos + i ] + '> '
                tags[  mask_pos + i  ] = 'X'

            else:
#        target_sequence[ snp_pos - pos + 1  ] = ' {' + target_sequence[ snp_pos - pos + 1  ] + '} '
                pass





    sequence =  "".join( sequence )
    sequence = re.sub(' ', '', sequence)



    tagged_string = "".join( tags )


    return ( sequence, tagged_string  )



def check_if_primer_clash( mapped_list, start, end):
    
    verbose_print( "check_if_primers_clash", 4)
    for i in range (start, end):
        if (mapped_list[ i ]  != " "):
            return True 

    return False


def extract_passed_primer_seqs( primer3_results, passed_primers):
    
    verbose_print( "extract_passed_primer_seqs", 3)

    primer_seqs = []
    for primer in passed_primers:
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        primer_seqs.append( [name, primer3_results[ primer ]] )


    return primer_seqs

def make_primer_mapped_strings( target_sequence, passed_primer_seqs):

    verbose_print( "make_primer_mapped_strings", 3)

    mappings = align_primers_to_seq(target_sequence, passed_primer_seqs )
    mapped_strings = []
    mapped_strings.append( [" "]*len( target_sequence ) )

    mapped_colours = []
    mapped_colours.append( [ -1 ]*len( target_sequence ) )

    for mapping in mappings:
        ( name, primer, pos, strand ) = mapping

        primer_nr = int(re.sub(r'.*_(\d)',r'\1' , name))



        mapping_index = 0

        arrows = (len(primer)-len(name)-2)/2
        arrow_type = ">"

        if ( strand == 1 ):
            primer = revDNA( primer )
            arrow_type = "<"

        tag = arrow_type*arrows + " " + name + " " + arrow_type*arrows
        if ( len(tag) < len(primer)):
            tag += arrow_type
            
        primer = tag 

        primer = list( primer )

        scan_index = 0
        while( 1 ):
            if (check_if_primer_clash ( mapped_strings[ scan_index ], pos, pos + len( primer ))):


#                print "increasing scan_index by one " + str( scan_index )
                import time
                time.sleep( 1 )
                scan_index += 1

                if (len( mapped_strings ) >= scan_index ):
                    mapped_strings.append([" "]*len( target_sequence ))
                    mapped_colours.append([ -1 ]*len( target_sequence ))

            else:
                break


        mapping_index = scan_index 

        for i in range(0, len( primer)):
            mapped_strings[ mapping_index ] [ pos + i ] = primer[ i ]
            mapped_colours[ mapping_index ] [ pos + i ] = primer_nr



    for i in range(0, len(mapped_strings)):
        mapped_strings[ i ] = "".join( mapped_strings[ i ])

    return mapped_strings, mapped_colours


def pretty_print_mappings( target_sequence, tagged_string, primer_strings, base1):

    lines = []
    
    for i in range(0, len(target_sequence), 80):


        lines.append("%-9d  %s" % ( base1+ i, target_sequence[i: i+80]))
        lines.append("           " + tagged_string[i: i+80])

        for primer_string in ( primer_strings ):

            line =  "           " + primer_string[i: i+80]
            if ( re.match(r'^ *$', line)):
                continue
            lines.append( line )

        lines.append( "" )

    return lines

    lines.append( "Map keys::" )
    lines.append( "XXXXXX excluded region" )
    lines.append( "****** target" )
    lines.append( ">>>>>> left primer" )
    lines.append( "<<<<<< right primer" )

    lines.append( "\n" )
    lines.append( "\n" )

    return lines


def pretty_pdf_mappings(top_offset,  target_sequence, tagged_string, primer_strings, primer_colours, base1):

    lines = []
    
    for i in range(0, len(target_sequence), 80):


#        c.drawString(40, top_offset, "%-9d  %s" % ( base1+ i, target_sequence[i: i+80]))
#        top_offset -= 8
#        c.drawString(40, top_offset, "           " + tagged_string[i: i+80])
#        top_offset -= 8

        p_line = "%-9d  %s" % ( base1+ i, target_sequence[i: i+80])
        m_line = "           " + tagged_string[i: i+80]


        x_offset = 40
        for k in range(0, len(p_line)):

            if (m_line[ k ] == "X"):
                c.setFillColorRGB(255,0,0)
            elif (m_line[ k ] == "*"):
                c.setFillColorRGB(0,190,0)

            c.drawString(x_offset , top_offset, p_line[k])
            x_offset += stringWidth(" ", 'mono', 8)
            c.setFillColorRGB(0,0,0)


        top_offset -= 8

        m_line = re.sub(r'X', r' ', m_line)

        if (re.search(r'\*', m_line)):
            c.setFillColorRGB(0,190,0)
            c.drawString(40 , top_offset, m_line)
            c.setFillColorRGB(0,0,0)

            top_offset -= 8


#        if (m_line.search(r'\*', name) :


        for j in range(0, len(primer_strings)):
            primer_string = primer_strings[ j ]
            primer_colour = primer_colours[ j ]

            line = primer_string[i: i+80]
            if ( re.match(r'^ *$', line)):
                continue

            x_offset = 40 + stringWidth(" ", 'mono', 8)*11

            for k in range(i, i+80):
                if ( k > len(target_sequence) - 1):
                    break

                if ( primer_colour[k] >= 0 ):

                    c.setFillColorRGB(colours[ primer_colour[k] ][0], 
                                      colours[ primer_colour[k] ][1],
                                      colours[ primer_colour[k] ][2])

                c.drawString(x_offset , top_offset, primer_string[k])
                x_offset += stringWidth(" ", 'mono', 8)
                c.setFillColorRGB(0,0,0)


            top_offset -= 8

        top_offset -= 8

#            lines.append( line )

#        lines.append( "" )


    return top_offset


#        c.drawString(40, top_offset, "Map keys::" )
#        top_offset -= 8
#        c.drawString(40, top_offset, "****** target" )
#        top_offset -= 8
#        c.drawString(40, top_offset, ">>>>>> left primer" )
#        top_offset -= 8
#        c.drawString(40, top_offset, "<<<<<< right primer" )
#        top_offset -= 8



def pretty_print_primer_data(primer3_results, passed_primers ):

    lines = []


    verbose_print( "extract_passed_primer_seqs", 3)

    lines.append( "\n" )
    lines.append( "\n" )
    lines.append( "\n" )

#    pp.pprint( primer3_results )

    lines.append( "_-=-"*15 +"_" )

    lines.append( " Primer design report for chr: %s pos: %d" % (chr, pos))
    lines.append( "_-=-"*15 +"_")
    lines.append( "\n")

#    lines.append( "\t".join(['ID', '%GC', 'TM', 'Sequence']))
    lines.append( "ID         %GC    TM     Primer sequence           Mapping(s)    ")
    lines.append( "_-=-"*15 +"_")

    primer_seqs = []
    for primer in sorted(passed_primers):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)


#        lines.append( "\t".join([name, 
#                         primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
#                         primer3_results[ "PRIMER_" + name + "_TM"],
#                         primer3_results[ "PRIMER_" + name + "_SEQUENCE"], 
#                         passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))


        lines.append( "%-10s %.2f  %.2f  %-25s %s" % (name, 
                                                     float(primer3_results[ "PRIMER_" + name + "_GC_PERCENT"]), 
                                                     float(primer3_results[ "PRIMER_" + name + "_TM"]),
                                                     primer3_results[ "PRIMER_" + name + "_SEQUENCE"], 
                                                     passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]))


    lines.append( "" )
    lines.append( "-="*46 )
    lines.append( "" )

    return lines



def pretty_pdf_primer_data(c, y_offset, primer3_results, passed_primers ):

#    c.drawString(40 , y_offset, "_-=-"*15 +"_" )
    c.line(40,y_offset,width - 40  ,y_offset+2)
    y_offset -= 8

    c.drawString(40 , y_offset, "Primer design report for chr: %s pos: %d" % (chr, pos))
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 16

    c.drawString(40 , y_offset, "ID         %GC    TM     Primer sequence           Mapping(s)    ")
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8




    primer_seqs = []
    for primer in sorted(passed_primers):
        if ( primer == 'FULLSEQ'):
            continue


        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        if (name == "RIGHT_0"):
            y_offset -= 8


#        lines.append( "\t".join([name, 
#                         primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
#                         primer3_results[ "PRIMER_" + name + "_TM"],
#                         primer3_results[ "PRIMER_" + name + "_SEQUENCE"], 
#                         passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))


        c.drawString(40 , y_offset, "%-10s %.2f  %.2f  %-25s %s" % ("", 
                                                     float(primer3_results[ "PRIMER_" + name + "_GC_PERCENT"]), 
                                                     float(primer3_results[ "PRIMER_" + name + "_TM"]),
                                                     primer3_results[ "PRIMER_" + name + "_SEQUENCE"], 
                                                     passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]))

        primer_nr = re.sub(r'.*_(\d)',r'\1' , name)

#        pp.pprint( primer_nr )
#        print pp.pprint(colours[ int( primer_nr ) ])
        

        c.setFillColorRGB(colours[ int( primer_nr ) ][0], 
                          colours[ int( primer_nr ) ][1],
                          colours[ int( primer_nr ) ][2])

        c.drawString(40 , y_offset, name)
        c.setFillColorRGB(0,0,0)
        y_offset -= 8


            



    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8
    y_offset -= 8


    return y_offset


def method_blurb():

    lines = []

    lines.append('primer-designer version: ' + VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

    return lines


def pretty_pdf_method(top_offset):

    lines = method_blurb()

    top_offset = 80

    for line in lines:
        c.drawString(40, top_offset, line)
        top_offset -= 8

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


args = get_and_parse_arguments()

chr = args.chr
pos = int(args.pos)
FLANK = args.flank or FLANK
FLANK = int( FLANK )

region_id = "%s:%d" % ( chr, pos)



target_sequence                   =   fetch_region( chr, pos - FLANK, pos + FLANK     )
(tagged_sequence, tagged_string)  =   markup_sequence( chr, pos, FLANK, target_sequence  )


primer3_results  = run_primer3( region_id , tagged_sequence, "%s_%d.primer3" % ( chr, pos))
passed_primers   = check_primers( region_id, target_sequence, primer3_results)
#pp.pprint( passed_primers )
passed_primer_seqs = extract_passed_primer_seqs( primer3_results, passed_primers )


(mapped_primer_strings, mapped_primer_colours) = make_primer_mapped_strings( target_sequence, passed_primer_seqs)



#lines  =  pretty_print_primer_data(primer3_results, passed_primers )

lines = pretty_print_mappings( target_sequence, tagged_string, mapped_primer_strings, pos - FLANK)



if (  args.text_output ):
    print "\n".join(lines)
    exit()

else:

#    lines += method_blurb()


#    print "\n".join(lines)

    filename = "%s_%d" % (chr, pos)

    if ( args.output ):
        filename = args.output + "_" + filename

    if (not re.match(".pdf", filename )):
        filename += ".pdf"

#    print "Going with filename '%s'" % filename


    c = canvas.Canvas( filename , pagesize=A4)
    
    width, height = A4 #keep for later
    print width

    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfbase.pdfmetrics import stringWidth 

    font = TTFont('mono', '/usr/share/fonts/truetype/ttf-liberation/LiberationMono-Regular.ttf')
    pdfmetrics.registerFont( font )
    c.setFont('mono', 7)
    top_offset = pretty_pdf_primer_data(c, height - 30, primer3_results, passed_primers )

    c.setFont('mono', 8)
    pretty_pdf_mappings(top_offset,  target_sequence, tagged_string, mapped_primer_strings, mapped_primer_colours, pos - FLANK)

    pretty_pdf_method(top_offset)



    # for line in lines:
    #     top_offset -= 8
    #     if (line == "\n"):
    #         continue

    #     if ( re.search(r'X', line) or re.search(r'\*', line)):
    #         x_offset = 40
    #         for character in line:
    #             if (character == 'X'):
    #                 c.setFillColorRGB(255,0,0)
    #                 c.drawString(x_offset , top_offset, character)
    #                 c.setFillColorRGB(0,0,0)
    #             elif  (character == '*'):
    #                 c.setFillColorRGB(0,255,0)
    #                 c.drawString(x_offset , top_offset, character)
    #                 c.setFillColorRGB(0,0,0)
    #             else:
    #                 c.drawString(x_offset , top_offset, character)

    #             x_offset += stringWidth(character, 'mono', 7)


    #     else:
    #         c.drawString(40 , top_offset, line)

    c.showPage()
    c.save()





