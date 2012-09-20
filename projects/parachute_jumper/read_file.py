#!/usr/bin/env python
import sys


def read_file(filenames):
	'''Method to read files and store variables them as dictionaries'''
	
	# Making sure the list is not empty
	if len(filenames) < 1:
		print 'Please make sure the list is not empty'
		sys.exit(1)
	
	dic={}
	for filename in filenames:
		try:
			f = open(filename,'r')
			for line in f:
				pair = line.split()
				dic[pair[0]] = eval(pair[1])
			f.close()
		except:
			print 'Something failed'
			sys.exit(1)

	return dic

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please include filenames at command line'
	else:
		mydict = read_file(sys.argv[1:])
	print mydict  


