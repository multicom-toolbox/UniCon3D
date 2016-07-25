############################################################################
# Author: Jie Hou, Debswapna Bhattacharya
# Date: 2016/07/25
############################################################################

import sys
import os
import optparse    # for option sorting
import string, re

parser = optparse.OptionParser()
parser.add_option('--fasta', dest='fasta',
	default = '',    # default empty!
	help = 'fasta file containing the target sequence')


(options,args) = parser.parse_args()

fasta = options.fasta

#open the fasta file
f = open(fasta, 'r')
text = f.readlines()    # read the text
f.close()    # close it
# removing the trailing "\n" and any header lines
text = [line.strip() for line in text if not '>' in line]
text = ''.join( text )    # combine into a single sequence
seqlen = len(text)

# Creates a list containing lists initialized to 0
CM = [[0 for x in range(seqlen)] for x in range(seqlen)] 

#now get ready to print it
for i in xrange(seqlen):
	for j in xrange(seqlen):
		sys.stdout.write(str(CM[i][j]) + ' ')
	sys.stdout.write('\n')
