#!/usr/bin/python 
# 
# 
# 
# 
# Kim Brugger (21 Oct 2015), contact: kbr@brugger.dk

import sys
import os
import pprint
pp = pprint.PrettyPrinter(indent=4)

import re

FLANK              = 500
NR_PRIMERS         = 4
ALLOWED_MISMATCHES = 4
MAX_MAPPINGS       = 5
MAX_PRODUCT_SIZE   = 800
MIN_PRODUCT_SIZE   = 120
 

def check_primers( smalt_results ):


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
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosomes: " + ",".join(set ( res[ primer ][ 'CHR' ] ))


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes



#    pp.pprint( smalt_report)
#    pp.pprint( res )


    return res




def digital_PCR( primer_mappings ):

#    pp.pprint( primer_mappings )


    primer_names = sorted(primer_mappings.keys())
    nr_primer_names = len( primer_names )

    mappings = {}

    products = {}
    

    for i in range(0, nr_primer_names):
        primer1 = primer_names[ i ] 

        if ( primer1 == 'FULLSEQ'):
            continue

        if ( not re.search(r'LEFT', primer1 )):
            continue


        mappings[ primer1 ] = {}
        products[ primer1 ] = {}

        for j in range(0, nr_primer_names):
            primer2 = primer_names[ j ] 

            if ( primer2 == 'FULLSEQ'):
                continue

            if ( not re.search(r'RIGHT', primer2 )):
                continue


            mappings[ primer1 ][ primer2 ] = []
            products[ primer1 ][ primer2 ] = []


#            print " -- %s vs %s"  % (primer1, primer2)


            for chr_index1 in range(0, len(primer_mappings[ primer1 ][ 'CHR' ])):
                for chr_index2 in range(0, len(primer_mappings[ primer2 ][ 'CHR' ])):

                    chr1 = primer_mappings[ primer1 ][ 'CHR' ][ chr_index1 ]
                    chr2 = primer_mappings[ primer2 ][ 'CHR' ][ chr_index2 ]

                    if ( chr1 != chr2 ):
                        continue


                    pos1 = int( primer_mappings[ primer1 ][ 'POS' ][ chr_index1 ] )
                    pos2 = int( primer_mappings[ primer2 ][ 'POS' ][ chr_index2 ] )

                    product_size = ( pos2 - pos1 )

#                    if ( product_size > MAX_PRODUCT_SIZE ):
#                        continue


#                    if ( product_size < 0 or 
#                         product_size > MAX_PRODUCT_SIZE ):
#                        continue

#                    if ( product_size < MIN_PRODUCT_SIZE ):
#                        continue

#                    print " %s:%d vs %s:%d ==> %d" % ( chr1, pos1, chr2, pos2, product_size )

                    mappings[ primer1 ][ primer2 ].append( product_size )
                    products[ primer1 ][ primer2 ].append( {'chr' : chr1, 'start_pos': pos1, 'end_pos': pos2, 'size': product_size} )

#        print "\n"



#    pp.pprint( products )
    return products
    exit()

    longest_product = 0
    longest_product_primer_pairs = ()
    for primer1 in mappings.keys():
        for primer2 in mappings[ primer1 ].keys():
        
            if ( len( mappings[ primer1 ][ primer2 ]) == 0 ):
                print "No usable pcr product from %s and %s" % ( primer1, primer2 )
                continue
            elif ( len( mappings[ primer1 ][ primer2 ]) > 1 ):
                print "multiple pcr products from %s and %s" % ( primer1, primer2 )
                continue

            print "%s + %s > %s bp " % (primer1, primer2,  mappings[ primer1 ][ primer2 ][ 0 ])

            if ( mappings[ primer1 ][ primer2 ][0] > longest_product ):
                longest_product = mappings[ primer1 ][ primer2 ][ 0 ]
                longest_product_primer_pairs = ( primer1, primer2 )
#                print "%s > %s (post)" % (longest_product,  mappings[ primer1 ][ primer2 ][ 0 ])



    print "\n\nLongest product (%d bp) comes from the %s and %s primer pair" % (longest_product, 
                                                                            longest_product_primer_pairs[0], 
                                                                            longest_product_primer_pairs[1])

                
def check_PCR_products(products, target_chr, target_start, target_end):


#    pp.pprint( products )

    failed_primer_pairs = []

    for primer1 in products:
        for primer2 in products[ primer1 ]:

            good_products = []
            bad_products  = []
            for product in products[ primer1 ][ primer2 ]:
#                pp.pprint( product )


                if ( product['end_pos'] - product['start_pos'] + 1 < MIN_PRODUCT_SIZE or
                     product['end_pos'] - product['start_pos'] + 1 > MAX_PRODUCT_SIZE):

#                    print "Wrong product size: %d " % ( product['end_pos'] - product['start_pos'] + 1 )
                    pass
#                    product = {}
                elif ( target_chr != product[ 'chr' ] ):
                    print " mis priming on diff chromosome " 
#                    products[ primer1 ][ primer2 ] = []
                    bad_products.append( product )

                elif( target_start < product['start_pos'] or target_end > product['end_pos'] ):
                    print "%s > %s or %s < %s " % ( target_start, product['start_pos'],  target_end, product['end_pos'] )
                    print "wrong region ( %s & %s )"  % ( primer1, primer2)
                    bad_products.append( product )
#                    products[ primer1 ][ primer2 ] = []
                else:
                    good_products.append( product )


            products[ primer1 ][ primer2 ] = { 'good': good_products, 'bad' : bad_products }



#    pp.pprint( products )

    return products

smalt_file = '8:96259936.smalt'

def pick_best_primers( primer_data, chromo, start_pos, end_pos ):
    
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
            print "Non unique mapping ( %s )" % primer
            continue


        if ( primer_data[ primer ][ 'CHR' ][ 0 ] != chromo ):
            print "Unique mapping to different chromosome (%s). Should never happen! " % primer
            continue

        if  (primer.find( 'LEFT' ) >= 0):
            
            primer_dist = start_pos - int (primer_data[ primer ][ 'POS' ][ 0 ]) + 1

            if ( primer_dist < 0 ):
                print "Primer %s downstream of region ! ( %d [%d, %d])" % (primer, primer_dist, start_pos, int (primer_data[ primer ][ 'POS' ][ 0 ]))
                continue

            if ( dist_fwd is None or primer_dist < dist_fwd):
                dist_fwd    = primer_dist
                closest_fwd = primer
                continue

        
        elif( primer.find( 'RIGHT' ) >= 0):

            primer_dist =  int (primer_data[ primer ][ 'POS' ][ 0 ]) - end_pos + 1

            if ( primer_dist < 0 ):
                print "Primer %s uptream of region ! (%d)" % (primer, primer_dist)
                continue

            if ( dist_rev is None or primer_dist < dist_rev ):
                dist_rev    = primer_dist
                closest_rev = primer
                continue

    

    return closest_fwd, closest_rev



if ( sys.argv >= 1 ):
    smalt_file = sys.argv[1]


region = smalt_file.rstrip(".smalt")
(chromo, pos) = region.split(":")
(start_pos, end_pos) = map(int, pos.split("-"))


primer_data = check_primers( smalt_file )

#pp.pprint( primer_data )


pcr_products = digital_PCR( primer_data )
pcr_products = check_PCR_products( pcr_products, chromo, start_pos, end_pos )


fwd_primer, rev_primer = pick_best_primers(primer_data, chromo, start_pos, end_pos)


print " Picked Primer Pair ( %s, %s )" % ( fwd_primer, rev_primer)


print "SMALT FILE :: %s " % smalt_file
