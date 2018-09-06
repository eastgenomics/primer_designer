#----------------------------------------------------------------
# Script for automating Sanger primer design. 
#----------------------------------------------------------------
# Nichols S. Gleadall & Kim Brugger
# 08/09/2015
# nick.gleadall@googlemail.com


#----------------------------------------------------------------
# Import modules
#----------------------------------------------------------------
import subprocess
import primer_3 as p3
import pipeliners as pipe
import pprint as pp


#----------------------------------------------------------------
# I/O (currently set example)
#----------------------------------------------------------------
sample_name = "Seq"
#sequence = "ATTTCTCATTCTAGCTTGCTCCCGCTCTGCAGTCTGTCCACAATCTTTCTCAAAGGGGGCAGGACAGGGAGCACCAATGTAAAATATCTTCCCTTTATTGCTCTTCTGTCTGTAACATACCCCAGCATAAGAAATATGCTATAAACATGAAGGAGGGCATGGGTGCAAAAAGGATGTGTTCATAAGTGCTCCAGCTCATAAGCCCAGTTCAAAATTGGGTCTAGCTAGACATAAGATGAACCAAATATCTCACCATTGATCAAGAAGAAGAAGGCTCCAGTCTCATTTGCTACAGCTCGAGCAATCAGGGTCTTTCCTGTTCCAGGAGGTCCGTAAAGCAGGATTCCTCTAGGAGGCTGTGGACCAATGCAGAGTGTGATTAGGTTCCCACCTTCTCCCAAACCCTAAAATGGCTTGACTAGCGCTCCAGAGAGGGTGAGAATCAGGCTTCAACCCTGACAATATTAAGAATAAGCCCTCCGCAGTAAATGCTAAGAAGAC"
sequence = "ACGAAGAGAAAAGGCCTCTCGCAGTCTGGAAATGGACACGTCTTTTTTTTCTTCTCCAGTCCTTAATGAAGGGGATTGATCCTGCTTTTCTACCATGGGCTTTTCCAAATCCGCTGCATGCATTTTTATTAAGTTACCTAAGCAAACGTGGACGGAGAAGAGGGTCAGGGACTATCCTGAAATGGTGAGAGGACGTGCTTATGTGAACAGATACTTCACAAAAGAGGAGATCCACATGCTAATTACACAGATGAACACAGTTCAATGTTCAAAATAAAACTATAATATGGGCCAGGTGTGGTGGCTTACGCCTGTTATCCCAGCACTTTAGGAGGCCAAGGCAGGGGGATCACATGAGGCTAGGAGTTCAGGACTGGTCTGGACAACATGGTAAAACCCTGTTTCTACTAAAAATACAAAAATTAGCCGGGTGTGGTGGCATATCTGTCATCCCAGCTACTTGGGAGGCTGAGGCACGAGAATCCCTTTAGC"
#----------------------------------------------------------------
# Define Functions
#----------------------------------------------------------------

# This captures the $_ output of the primer3 engine and returns it as a dictionary.
def run_p3( cmd ):

  process = subprocess.Popen(proc1_cmd,
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


#----------------------------------------------------------------
# System call / Subprocess commands
#----------------------------------------------------------------

# Create primer3 input document using primer_3 module
p3.make_input_file(sample_name,sequence)
	
# Run the primer3 engine using the function above and get back the results in a dict. 
# proc1_cmd runs primer3 with -strict_tags set and the appropriate sample.prime3 file.
proc1_cmd = "/software/bin/primer3_core -strict_tags < " + sample_name + ".prime3"
p3_dict = run_p3( proc1_cmd )
#pp.pprint(p3_dict)

# print results to simple fasta

with open( sample_name + "_p3seq.fasta", 'w+') as primerfasta:

  primerfasta.write(">" + sample_name + "FULLSEQ\n" + p3_dict[ 'SEQUENCE_TEMPLATE' ] + "\n")

  primerfasta.write(">" + sample_name + "_1f\n" + p3_dict[ 'PRIMER_LEFT_0_SEQUENCE' ] + "\n")
  primerfasta.write(">" + sample_name + "_1r\n" + p3_dict[ 'PRIMER_RIGHT_0_SEQUENCE' ] + "\n")

  primerfasta.write(">" + sample_name + "_2f\n" + p3_dict[ 'PRIMER_LEFT_1_SEQUENCE' ] + "\n")
  primerfasta.write(">" + sample_name + "_2r\n" + p3_dict[ 'PRIMER_RIGHT_1_SEQUENCE' ] + "\n")

  primerfasta.write(">" + sample_name + "_3f\n" + p3_dict[ 'PRIMER_LEFT_2_SEQUENCE' ] + "\n")
  primerfasta.write(">" + sample_name + "_3r\n" + p3_dict[ 'PRIMER_RIGHT_2_SEQUENCE' ] + "\n")

  primerfasta.write(">" + sample_name + "_4f\n" + p3_dict[ 'PRIMER_LEFT_3_SEQUENCE' ] + "\n")
  primerfasta.write(">" + sample_name + "_4r\n" + p3_dict[ 'PRIMER_RIGHT_3_SEQUENCE' ] + "\n")

  primerfasta.write(">" + sample_name + "_5f\n" + p3_dict[ 'PRIMER_LEFT_4_SEQUENCE' ] + "\n")
  primerfasta.write(">" + sample_name + "_5r\n" + p3_dict[ 'PRIMER_RIGHT_4_SEQUENCE' ] + "\n")

primerfasta.close()


# Re define primer fasta file name for smalt use.
primerfasta = sample_name + "_p3seq.fasta"
# Run SMALT using g1kv37 as ref and primer fasta's as query seq's - USED HASH INDEX k13 and S5.
proc2_cmd = "/software/bin/smalt map -d -1 /refs/human_1kg/human_g1k_v37 " + primerfasta + "> " + sample_name + ".smalt"
subprocess.call(proc2_cmd, shell=True)

smalt_output = sample_name + ".smalt"
# Check the SMALT results - ensure that primers match well and are within the same region as origional sequence. 

id_word = "FULLSEQ"
variant_region = {}
smalt_report = {}
query_region = []

with open(smalt_output, 'rU') as primer_check:
  for line in primer_check:
    
    if not line.startswith("cigar"):
      continue
       
    field = line.split(" ")
    variant_region[ 'name'       ] = field[  1  ] #mapping_score
    variant_region[ 'chromosome' ] = field[  5  ] #mapping_score
    variant_region[ 'startpos'   ] = field[  6  ] #mapping_score
    variant_region[ 'endpos'     ] = field[  7  ] #mapping_score
    variant_region[ 'len'        ] = field[  9  ] #mapping_score
    variant_region[ 'lenmatched' ] = field[  11 ] #mapping_score

    smalt_report.append(variant_region)
    


print smalt_report


    




    