#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (11 Jul 2016), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

import os
import sys
import re
import string
import random
#import pipeliners as pipe

import subprocess
def system_call( cmd ):

    try:
        subprocess.check_call(cmd, shell=True)

    except subprocess.CalledProcessError as scall:
        print "System call '%s' failed. Exit code: %s Error output: %s" %( cmd,  scall.returncode, scall.output)
        exit()



# Define system call

if ( len(sys.argv ) == 1 ):
	print "USAGE: bulk_design.py input-file [working dir]"
	print "input-file one tab seperated output-name, region and reference (grch37 or grch38) of interest pr line"
	exit()

infile = sys.argv[1]
if ( infile.find(".txt") == -1):
	print "Infile '%s' does not look like a txt file" % infile
	exit()


if ( len(sys.argv ) == 3 ):
    os.chdir( sys.argv[ 2 ])


outfile = re.sub(r'.txt', '.zip', infile)

# Run designer function 


def random_string(length=10):
    """ Creates a random string 

    Args:
      length (int): length of string to generate, default 10

    returns:
      string

    raises:
      None
    """

    random_string = ""
    # Choose from lowercase, uppercase,  and digits
    alphabet = string.ascii_lowercase + string.ascii_uppercase + string.digits
    for n in range(0, length):
        random_string += random.choice( alphabet )

    return random_string

if __name__ == '__main__':


    working_dir = random_string()
    os.mkdir( working_dir )
    os.chdir( working_dir )


    # Main Loop
    with open( "../{}".format( infile ) , "rU" ) as plist:
	
	for line in plist:

            print line


            line = line.strip("\n")
            if ( line == '' ):
                continue 

            (testname, region, reference) = re.split(r"[\t ]+", line)

            reference = reference.lower()

            (chrom, npos) = region.split(":")
            cmd = "/mnt/storage/apps/software/primer_designer/1.1/primer_designer_region.py -c {} -p {} -o {} --{}".format(chrom, npos, testname, reference)
            print cmd 
            system_call(cmd )


    os.chdir( ".." )
    system_call( "zip -j {} {}/*pdf".format( outfile, working_dir ))

    print "SUCCESS\nOutput file: %s" % outfile
    
