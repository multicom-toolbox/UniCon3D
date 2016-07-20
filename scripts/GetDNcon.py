#!usr/bin/env python


############################################################################
# Author: Jie Hou, Debswapna Bhattacharya
# Date: 2016/07/20
############################################################################


import sys
import os
import operator 
import optparse    # for option sorting 



dncon_exe = '/home/jh7x3/tools/dncon1.0/run_dncon.pl'
dncon_submit = '/home/jh7x3/tools/dncon1.0/local_submit_job.cgi'

class DNcon():
	"""
	The DNcon class
	"""
	def __init__(self, fasta_file, contact_range,CPU):
		#constructor, initialzing a couple values
		self.fasta_file = fasta_file		
		self.contact_range = contact_range
                self.CPU = CPU
                self.curr_dir = os.getcwd()
		self.dncon_dir = self.curr_dir + '/DNcon_contact'
		if not os.path.exists(self.dncon_dir):
			os.makedirs(self.dncon_dir)
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
		os.chdir(self.dncon_dir)

                FASTA_dir = self.dncon_dir + '/FASTA'
                if not os.path.exists(FASTA_dir):
                        os.makedirs(FASTA_dir)
                my_cmd = 'cp ' + self.fasta_file + ' ' + FASTA_dir + '/'
                retcode = os.system(my_cmd)
                if retcode != 0:
                        print 'copy fasta: Failed to copy fasta'
                        sys.exit(1)

                Output_dir = self.dncon_dir + '/output'
                if not os.path.exists(Output_dir):
                        os.makedirs(Output_dir)
                log_file = os.path.join(self.log_dir,'DNcon.log')
		my_cmd = 'perl ' + dncon_exe + ' ' + dncon_submit  + ' ' + FASTA_dir + ' ' + Output_dir + ' ' + str(self.contact_range) + ' ' + str(self.CPU) + ' &> ' + log_file
		retcode = os.system(my_cmd)
		if retcode != 0:
			print 'DNcon: Failed to predict'
			sys.exit(1)
		os.chdir(self.curr_dir)

parser = optparse.OptionParser()
parser.add_option('--fasta', dest='fasta',
        default = '',    # default empty!
        help = 'fasta file containing the target sequence')

parser.add_option('--range', dest='range',
        default = '',    # default empty!
        help = 'short-range_contact-num(l_5,l_10,l,2l,5l,10l)')

parser.add_option('--device', dest='device',
        default = '',    # default empty!
        help = 'Device: CPU/GPU')

(options,args) = parser.parse_args()

fasta_file = options.fasta
contact_range = options.range
device = options.device


print 'predicting contact prediction from DNcon ...'
step = DNcon(fasta_file,contact_range,device)
step.apply()
print '              done'
print ''


