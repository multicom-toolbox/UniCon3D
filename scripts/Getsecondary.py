#!usr/bin/env python

############################################################################
# Author: Jie Hou, Debswapna Bhattacharya
# Date: 2016/07/20
############################################################################


import sys
import os
import optparse    # for option sorting

scratch_exe = '/home/jh7x3/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh'

class secondary():
	"""
	The secondary class
	"""
	def __init__(self, fasta_file):
		#constructor, initialzing a couple values
		self.fasta_file = fasta_file		
		self.target = 'test'
		self.cpu = 15
		self.curr_dir = os.getcwd()
		self.secondary_dir = self.curr_dir + '/secondary'
		if not os.path.exists(self.secondary_dir):
			os.makedirs(self.secondary_dir)
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
		os.chdir(self.secondary_dir)
		log_file = os.path.join(self.log_dir,'scratch.log')
		my_cmd = scratch_exe + ' ' + self.fasta_file + ' ' + self.target + ' ' + str(self.cpu) + ' &> ' + log_file
		retcode = os.system(my_cmd)
		if retcode != 0:
			print 'secondary: Failed to predict'
			sys.exit(1)
		os.chdir(self.curr_dir)

parser = optparse.OptionParser()
parser.add_option('--fasta', dest='fasta',
        default = '',    # default empty!
        help = 'fasta file containing the target sequence')
parser.add_option('--id', dest='id',
        default = '',    # default empty!
        help = 'target id')

(options,args) = parser.parse_args()

fasta_file = options.fasta
target = options.id

print 'predicting secondary structure...'
step2 = secondary(fasta_file)
step2.set_target(target)
step2.apply()
print '              done'
print ''


