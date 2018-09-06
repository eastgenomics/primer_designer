#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (15 Sep 2015), contact: kim@brugger.dk

import sys
sys.stdout.flush()
import os
import pprint
pp = pprint.PrettyPrinter(indent=4)
import argparse
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import stringWidth 

import shlex
import subprocess
import re

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
width, height = A4 #keep for later

TMP_FILES = []

colours = [[255,   0,   0], # red
           [  0, 255,   0], # green
           [  0,   0, 255], # blue
           [255,   0, 255], # Pink
           [0,   255, 255],
           [255, 255,   0]] # Yellow, crap!

# External tools we use 
SAMTOOLS   = '/mnt/storage/apps/software/samtools/1.5/bin/samtools ';
REFERENCE  = None
DBSNP      = None
TABIX      = '/mnt/storage/apps/software/samtools/1.5/bin/tabix '
SMALT      = '/mnt/storage/apps/software/smalt/0.7.6/bin/smalt '
PRIMER3    = '/mnt/storage/apps/software/primer3_core/2.3.7/src/primer3_core '

# default parameters, can be overridden by the user (perhaps)
TARGET_LEAD        = 50
FLANK              = 500
NR_PRIMERS         = 4
ALLOWED_MISMATCHES = 4
MAX_MAPPINGS       = 5

VERBOSE    =  2
VERSION    =  '1.1'

def verbose_print( msg, level ):
    if ( level <= VERBOSE ):
        print msg

def fetch_region( chr, start, end ):

    verbose_print( "fetch_region", 2)

    cmd = "%s faidx %s  %s:%d-%d " % ( SAMTOOLS, REFERENCE, chr, start, end )
#    print "cmd: {}".format( cmd )
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
    verbose_print( "fetch_known_SNPs", 2)

    cmd = "%s %s  %s:%d-%d " % ( TABIX, tabix_file, chr, start, end )

    verbose_print( cmd, 3)

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

#    pp.pprint( output );

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 



def write_primer3_file(seq_id, seq, primer3_file = ""):
    verbose_print( "write_primer3_file", 2)

    upstream_end      = 1
    downstream_start =  2
    for i in range(0, len(seq)):
            if seq[ i ] == '[':
                upstream_end = i
            elif seq[ i ] == ']':
                downstream_start = i


    template = '''SEQUENCE_ID={seq_id}
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
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/storage/apps/software/primer3_core/2.3.7/src/primer3_config/
='''

    details = {
        "seq_id" : seq_id,
        "seq"    : seq,
        "flank"  : upstream_end,
        "len"    : downstream_start - upstream_end + 1
    }

 
    with open( primer3_file, 'w+') as outfile:
        outfile.write(template.format(**details))

    outfile.close()


# This captures the $_ output of the primer3 engine and returns it as a dictionary.
def run_primer3( seq_id, seq, primer3_file = ""):
    verbose_print( "run_primer3", 2)

    if ( primer3_file == "" ):
        primer3_file = seq_id + ".primer3"


    primer3_file = re.sub(":", "_", primer3_file)


    TMP_FILES.append( primer3_file )

    write_primer3_file(seq_id, seq, primer3_file)


    cmd = PRIMER3 + " < " + primer3_file


    verbose_print( cmd, 3)
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
    verbose_print( "check_primers", 2)


    primers_file = region_id + "_p3seq.fasta"

    TMP_FILES.append( primers_file )

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
    TMP_FILES.append( smalt_results )

    ref = REFERENCE
    ref = re.sub(r'\..*', '', ref)

    cmd = "{} map  -d -1 {} {} > {}  ".format( SMALT, ref, primers_file, smalt_results)
    
    print cmd

#    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"
    #cmd = SMALT + " map  -d -1 /mnt/scratch/refs/homo_sapiens/GRCh37/Homo_sapiens_assembly37 " + primers_file + "> " + smalt_results + " 2> /dev/null"


    verbose_print( cmd, 4)
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
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosomes: " + ",".join(set(res[ primer ][ 'CHR' ]))


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes



#    pp.pprint( smalt_report)
#    pp.pprint( res )
#    exit()


    return res

def pick_best_primers( primer_data, chromo, start_pos, end_pos ):
    verbose_print("pick_best_primers", 2)

    # Finally we are getting to the crux of this whole ordeal: Pick the best primers.
    # I will be done using the following rules:
    # Unique MAPPING primers are best
    # primers closest to the region of interest
    # Primers generating odd products are eliminated.


    # First group the primers according to the region.

#    pp.pprint( primer_data )
#    pp.pprint( pcr_products )

    (closest_fwd, dist_fwd) = (None, None)
    (closest_rev, dist_rev) = (None, None)
    

    for primer in primer_data:

        if ( primer == 'FULLSEQ'):
            continue

        if ( primer_data[ primer][ 'MAPPING_SUMMARY' ] != 'unique mapping'):
            verbose_print( "Non unique mapping ( %s )" % primer, 5)
            continue

        if ( not primer_data[ primer ][ 'CHR' ] or primer_data[ primer ][ 'CHR' ][ 0 ] != chromo ):
            verbose_print( "No mapping or Unique mapping to different chromosome (%s). Should never happen! " % primer, 5)
            continue

        if  (primer.find( 'LEFT' ) >= 0):
            
            primer_dist = start_pos - int (primer_data[ primer ][ 'POS' ][ 0 ]) + 1

            if ( primer_dist < 0 ):
                verbose_print("Primer %s downstream of region ! ( %d [%d, %d])" % (primer, primer_dist, start_pos, int (primer_data[ primer ][ 'POS' ][ 0 ])), 5)
                continue

            if ( dist_fwd is None or primer_dist < dist_fwd):
                dist_fwd    = primer_dist
                closest_fwd = primer
                continue

        
        elif( primer.find( 'RIGHT' ) >= 0):

            primer_dist =  int (primer_data[ primer ][ 'POS' ][ 0 ]) - end_pos + 1

            if ( primer_dist < 0 ):
                verbose_print( "Primer %s uptream of region ! (%d)" % (primer, primer_dist), 5)
                continue

            if ( dist_rev is None or primer_dist < dist_rev ):
                dist_rev    = primer_dist
                closest_rev = primer
                continue

    

    return closest_fwd, closest_rev




def revDNA( string ):
    verbose_print( "revDNA", 4)
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
    verbose_print( "align_primers_to_seq", 2)

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


    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chr')
    
    # NICK ADDED RANGE
    
    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument('-p', '--pos')
    group.add_argument('-r', '--range', nargs = 2)

    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank')
    parser.add_argument('-t', '--text_output')

    parser.add_argument('--hg19', help='use hg19 as reference genome',action="store_true" )
    parser.add_argument('--grch37', help='use GRCh37 as reference genome',action="store_true")
    parser.add_argument('--grch38', help='use GRCh38 as reference genome',action="store_true")


    args = parser.parse_args()

    global REFERENCE
    global DBSNP


    if ( args.grch37 or args.hg19):
        REFERENCE = '/mnt/scratch/refs/human_1kg/human_g1k_v37.fasta'
	#REFERENCE = '/mnt/scratch/refs/homo_sapiens/GRCh37/Homo_sapiens_assembly37.fasta'
        DBSNP     = '/mnt/scratch/refs/homo_sapiens/GRCh37/annots-rsIDs-dbSNPv144.20150605.tab.gz'
    elif ( args.grch38 ):
        REFERENCE = '/mnt/scratch/refs/homo_sapiens/GRCh38/Homo_sapiens_assembly38.fasta'
        DBSNP     = '/mnt/scratch/refs/homo_sapiens/GRCh37/dbsnp_150_common.vcf.gz'
    else:
        print "Please select a reference genome to use"
        parser.parse_args('-h')

        
    


    return args



def markup_sequence( chr, startpos, endpos, flank, sequence):
    verbose_print( "tag_sequence", 2)

    sequence = list(sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    start = FLANK
    end = len(tags) - flank

    for x in range(0, len(tags)):
        if (x >= start and x <= end):
            tags[x] = '*'

        if x in range(start - TARGET_LEAD, start) or x in range(end, end + TARGET_LEAD):
            tags[x] = '-'
    #return tags
        

    verbose_print( "::::: %d - %d, %d" % (startpos, endpos, flank), 5)
    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ flank - TARGET_LEAD ] = ' [' + target_sequence [ flank - TARGET_LEAD ]
    sequence[(- flank + TARGET_LEAD) -1 ] = sequence [ (- flank + TARGET_LEAD) -1 ] + '] '

#    return sequence
#   exit ()

    dbSNPs = fetch_known_SNPs( DBSNP, chr, startpos - flank, endpos + flank )

    masked_positions = []

    for dbSNP in (dbSNPs):
        #print dbSNP
        #exit()
        if ( len( dbSNP ) < 6):
            continue
        
        #unpack dbSNP entry. 
        (snp_chr, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP

        snp_pos = int( snp_pos )
        if (snp_pos >= startpos - TARGET_LEAD and snp_pos <= endpos + TARGET_LEAD ):
            continue

        if ( common == '1'):
            #        pp.pprint( dbSNP )
            mask_pos = snp_pos - (startpos - flank) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
            
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


    #print tags
    #exit()
    tagged_string = "".join( tags )


    return ( sequence, tagged_string  )



def check_if_primer_clash( mapped_list, start, end):
    
    verbose_print( "check_if_primers_clash", 4)
    for i in range (start, end):
        if (mapped_list[ i ]  != " "):
            return True 

    return False


def extract_passed_primer_seqs( primer3_results, passed_primers):
    
    verbose_print( "extract_passed_primer_seqs", 2)

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

    verbose_print( "make_primer_mapped_strings", 2)

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
    verbose_print("pretty_pdf_mappings", 2)

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

        if (re.search(r'\*', m_line) or re.search(r'\-', m_line)):
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


    if ( startpos == endpos ):
        lines.append( " Primer design report for chr: %s position: %d" % (chr, startpos))
    else:
        lines.append( " Primer design report for chr: %s range: %d-%d" % (chr, startpos, endpos))

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



def pretty_pdf_primer_data(c, y_offset, primer3_results, passed_primers, fwd_primer=None, rev_primer=None ):
    verbose_print("pretty_pdf_primer_data", 2)

#    c.drawString(40 , y_offset, "_-=-"*15 +"_" )
    c.line(40,y_offset,width - 40  ,y_offset+2)
    y_offset -= 8


    if ( startpos == endpos ):
        c.drawString(40 , y_offset, "Primer design report for chr: %s position: %d" % (chr, startpos))
    else:
        c.drawString(40 , y_offset, "Primer design report for chr: %s range: %d-%d" % (chr, startpos, endpos))
            
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 16

    c.drawString(40 , y_offset, "ID         %GC    TM     Primer sequence           Best primer  Mapping(s)    ")
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

        picked_primer = ' '

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'


        c.drawString(40 , y_offset, "%-10s %.2f  %.2f  %-25s %s            %s" % ("", 
                                                     float(primer3_results[ "PRIMER_" + name + "_GC_PERCENT"]), 
                                                     float(primer3_results[ "PRIMER_" + name + "_TM"]),
                                                     primer3_results[ "PRIMER_" + name + "_SEQUENCE"], picked_primer,
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

def pretty_primer_data(outfile, primer3_results, passed_primers, fwd_primer=None, rev_primer=None ):
    verbose_print("pretty_primer_data", 2)


    fh = open( outfile, 'w')


    lines = []

    if ( startpos == endpos ):
        lines.append("Primer design report for chr: %s position: %d" % (chr, startpos))
    else:
        lines.append("Primer design report for chr: %s range: %d-%d" % (chr, startpos, endpos))

    lines.append("ID\t%GC\tTM\tPrimer sequence\tBest primer\tMapping(s)")

    for primer in sorted(passed_primers):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        picked_primer = ''

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'

        lines.append("\t".join([name, 
                              primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
                              primer3_results[ "PRIMER_" + name + "_TM"],
                              primer3_results[ "PRIMER_" + name + "_SEQUENCE"], picked_primer,
                              passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))




    lines.append("\n")
    fh.write("\n".join( lines ))
    fh.close()



def method_blurb(args):

    lines = []

    if args.grch37 or args.hg19:
        lines.append('primer-designer version: ' + VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')
    else:
        lines.append('primer-designer version: ' + VERSION + ' using dbSNP 150 for SNP checking, and human reference GRCh38.')
        

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')


    return lines


def pretty_pdf_method(top_offset, args):
    verbose_print("pretty_pdf_method", 3)

    lines = method_blurb(args)

    top_offset = 80

    for line in lines:
        c.drawString(40, top_offset, line)
        top_offset -= 8

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if ( __name__ == '__main__'):


    args = get_and_parse_arguments()

    chr = args.chr


    if ( args.range ):
        startpos, endpos = [int(x) for x in args.range]
    else:
        startpos = int( args.pos )
        endpos   = int( args.pos )


    FLANK = args.flank or FLANK
    FLANK = int( FLANK )
    
    if ( startpos == endpos):
        region_id = "%s:%d" % ( chr, startpos)
    else:
        region_id = "%s:%d-%d" % ( chr, startpos, endpos)



    target_sequence                   =   fetch_region( chr, startpos - FLANK, endpos + FLANK     )
#    print "Target seq: {}".format( target_sequence )
    (tagged_sequence, tagged_string)  =   markup_sequence( chr, startpos, endpos, FLANK, target_sequence  )



    primer3_results  = run_primer3( region_id , tagged_sequence, "%s_%d.primer3" % ( chr, startpos))
    passed_primers   = check_primers( region_id, target_sequence, primer3_results)
    passed_primer_seqs = extract_passed_primer_seqs( primer3_results, passed_primers )

    fwd_primer, rev_primer = pick_best_primers(passed_primers, chr, startpos, endpos)



    (mapped_primer_strings, mapped_primer_colours) = make_primer_mapped_strings( target_sequence, passed_primer_seqs)



    lines = pretty_print_mappings( target_sequence, tagged_string, mapped_primer_strings, startpos - FLANK)


    if (startpos == endpos):
        filename = "%s_%d" % (chr, startpos)
    else:
        filename = "%s_%d_%d" % (chr, startpos, endpos)


    if ( args.output ):
        filename = args.output + "_" + filename
        
    filename = filename.rstrip(".pdf")
    filename = filename.rstrip(".txt")

    if (1 or   args.text_output ):
        pretty_primer_data("%s.txt" % filename, primer3_results, passed_primers, fwd_primer, rev_primer )





    if (  args.text_output ):
        print "\n".join(lines)
        exit()
    
    else:

        filename = "%s.pdf" % filename


        c = canvas.Canvas( filename , pagesize=A4)
    
        width, height = A4 #keep for later


        font = TTFont('mono', '/usr/share/fonts/liberation/LiberationMono-Regular.ttf')
        pdfmetrics.registerFont( font )
        c.setFont('mono', 7)
        top_offset = pretty_pdf_primer_data(c, height - 30, primer3_results, passed_primers, fwd_primer, rev_primer )

        c.setFont('mono', 8)
        pretty_pdf_mappings(top_offset,  target_sequence, tagged_string, mapped_primer_strings, mapped_primer_colours, startpos - FLANK)

        pretty_pdf_method(top_offset, args)

        c.showPage()
        c.save()


    for filename in     TMP_FILES:
        print "deleting tmp file: %s " % filename 
#        os.remove( filename )






