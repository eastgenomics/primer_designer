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


chr = sys.argv[1]
pos = int(sys.argv[2])

SAMTOOLS  = '/software/bin/samtools ';
REFERENCE = '/refs/human_1kg/human_g1k_v37.fasta '
TABIX     = '/software/bin/tabix '

FLANK     = 250

def fetch_region( chr, start, end ):

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

    cmd = "%s %s  %s:%d-%d " % ( TABIX, tabix_file, chr, start, end )

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

#    pp.pprint( output );

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 



clean_target_sequence =  fetch_region( chr, pos - FLANK, pos + FLANK)
target_sequence       =  list( clean_target_sequence )
clean_target_sequence =  list( clean_target_sequence )
#target_sequence = target_sequence.split()

clean_target_sequence[ 5 ] = 'N'

#print len(target_sequence)

# Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
target_sequence[ FLANK ] = ' [' + target_sequence [ FLANK ] + '] '

dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/test.tab.gz', chr, pos - FLANK, pos + FLANK )

masked_positions = []

for dbSNP in (dbSNPs):
    if ( len( dbSNP ) < 6):
        continue

    (snp_chr, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP

    snp_pos = int( snp_pos )
    if ( common == '1'):
#        pp.pprint( dbSNP )
        mask_pos = snp_pos - (pos - FLANK) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
        if ( mask_pos == FLANK ):
            continue
        
        # In the odd case we are looking at common deletion, mask the whole region. Normally this will just be one base
        for i in range(0, len(snp_ref)):

            # If already masked skip masking it again.
            if ( re.search('<', target_sequence[ mask_pos + i  ]) or 
                 re.search('\[', target_sequence[ mask_pos + i  ])):
                continue

            target_sequence[ mask_pos + i  ] = ' <' + target_sequence[ mask_pos + i ] + '> '
            masked_positions.append( mask_pos + i )

    else:
#        target_sequence[ snp_pos - pos + 1  ] = ' {' + target_sequence[ snp_pos - pos + 1  ] + '} '
        pass





target_sequence =  "".join( target_sequence )
target_sequence = re.sub(' ', '', target_sequence)

print target_sequence


def revDNA( string ):
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



def align_primers_to_seq( seq, *all_primers):

    primers = []
    rev_primers = []
    for primer in all_primers:
        primers.append(  primer )
        rev_primers.append( revDNA( primer ))

    pp.pprint( primers )

    mappings = []

    for i in range(0, len(seq)):
        for primer in primers:
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
                print primer + " Matches at pos " + str( i )

                mappings.append( [primer, i, 0] )

        for primer in rev_primers:
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
                print primer + " Matches at pos " + str( i ) + " minus! "

                mappings.append( [primer, i, 1] )


    return mappings

primer_1_fwd = "TGTTATGGGTGAAACTCTGGGA"
primer_1_rev = "AAACAAACACAACCTGGAAGAC"

primer_2_fwd = 'ACTGATGTCTTGACTCATGGGT'
primer_2_rev = 'ACCTTGTGATATGTTTGCAGACA'



mappings = align_primers_to_seq(clean_target_sequence, primer_1_fwd, primer_1_rev, primer_2_fwd, primer_2_rev)

mapping_string = [" "] * len( clean_target_sequence )
mapping_string[ FLANK ] = '*'

for mapping in mappings:
    ( primer, pos, strand ) = mapping


    if ( strand == 1 ):
        primer = revDNA( primer )
        primer = "<" + primer + "<"
    else:
        primer = ">" + primer +">"

    pos -= 1

    primer = list( primer )


    for i in range(0, len( primer)):
        mapping_string[ pos + i ] = primer[ i ]

for masked_position in masked_positions:
    mapping_string[ int( masked_position) ] = 'X'


                        
pp.pprint( masked_positions )

print "".join( clean_target_sequence )
print "".join( mapping_string )
