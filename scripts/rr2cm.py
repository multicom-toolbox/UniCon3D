import sys
import os
import copy
import gzip
import math
import types
import threading
import datetime
import time
from collections import defaultdict
import string, re
import itertools as it
import itertools
import random
import linecache
import optparse    # for option sorting
import math

parser = optparse.OptionParser()
parser.add_option('--fasta', dest='fasta',
	default = '',    # default empty!
	help = 'fasta file containing the target sequence')
parser.add_option('--rr', dest='rr',
	default = '',    # default empty!
	help = 'rr file in CASP format')
parser.add_option('--xl', dest='xl',
	default = 1000000.0,    # default 100000000.0
	help = 'fraction of top contacts to select as a function of sequence length')


(options,args) = parser.parse_args()

fasta = options.fasta
rr = options.rr
xl = float(options.xl)

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

maxnum = seqlen * xl

num = 0
#open the rr file
r = open(rr, 'r')
for line in r:
	temp = re.split(" +",line)
	if len(temp) == 5:
		num += 1
		if num > maxnum:
			continue
		i = int(temp[0].strip()) - 1
		j = int(temp[1].strip()) - 1
		p = float(temp[4].strip())
		CM[i][j] = p
		CM[j][i] = p
r.close()    # close it		

#now get ready to print it
for i in xrange(seqlen):
	for j in xrange(seqlen):
		sys.stdout.write(str(CM[i][j]) + ' ')
	sys.stdout.write('\n')
