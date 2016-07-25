############################################################################
# Author: Jie Hou, Debswapna Bhattacharya
# Date: 2016/07/25
############################################################################

import sys
import os
import optparse    # for option sorting
import string, re

parser = optparse.OptionParser()
parser.add_option('--rr', dest='contact',
	default = '',    # default empty!
	help = 'residue-residue contact file with three column format (id id weight)')


(options,args) = parser.parse_args()

fasta = options.fasta

#open the fasta file
f = open(fasta, 'r')
text = f.readlines()    # read the text
f.close()    # close it
# removing the trailing "\n" and any header lines
text = [line.strip() for line in text if not '>' in line]
for line in text:
    print line.split( )

#now get ready to print it
#for i in xrange(seqlen):
#	for j in xrange(seqlen):
#		sys.stdout.write(str(CM[i][j]) + ' ')
#	sys.stdout.write('\n')
