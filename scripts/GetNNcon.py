#!usr/bin/env python

############################################################################
# Author: Jie Hou, Debswapna Bhattacharya
# Date: 2016/07/20
############################################################################

import sys
import os
import operator 
import optparse    # for option sorting 


NNcon_exe = '/home/jh7x3/tools/nncon1.0/bin/predict_ss_sa_cm.sh '

class NNcon():
	"""
	The NNcon class
	"""
	def __init__(self, fasta_file):
		#constructor, initialzing a couple values
		self.fasta_file = fasta_file		
                self.curr_dir = os.getcwd()
		self.NNcon_dir = self.curr_dir + '/NNcon_contact'
		if not os.path.exists(self.NNcon_dir):
			os.makedirs(self.NNcon_dir)
		self.log_dir = self.curr_dir + '/log'
		if not os.path.exists(self.log_dir):
			os.makedirs(self.log_dir)
	def __str__(self):
		"""
		All objects MUST include this method. 
		"""
		return self.get_name()
	def get_name(self): 
		""" 
		Return the name of the class of the object instance.
		""" 
		return self.__class__.__name__
	def set_target(self,target):
		"""
		This method sets the target
		"""
		self.target = str(target)
	def apply(self):
		"""
		This method performs the job
		"""
		os.chdir(self.NNcon_dir)
                log_file = os.path.join(self.log_dir,'NNcon.log')
		my_cmd = 'perl ' + NNcon_exe + ' ' + self.fasta_file + ' ' + self.NNcon_dir + ' &> ' + log_file
		retcode = os.system(my_cmd)
		if retcode != 0:
			print 'NNcon: Failed to predict'
			sys.exit(1)
		os.chdir(self.curr_dir)

parser = optparse.OptionParser()
parser.add_option('--fasta', dest='fasta',
        default = '',    # default empty!
        help = 'fasta file containing the target sequence')

(options,args) = parser.parse_args()

fasta_file = options.fasta

print 'predicting contact map by NNcon ...'
step2 =  NNcon(fasta_file)
step2.apply()
print '              done'
print ''


