#!/usr/bin/env python
import sys


def read_file(infilename0, infilename1):
	'''Method to read 2 files and store them as dictionaries'''
	try:
	    infilename0 = sys.argv[1]; infilename1 = sys.argv[2]
	except:
	    "Usage:",sys.argv[0], "infile0 infile1";sys.exit(1)
	
	infile0 = open(infilename0, "r") # Open first file for reading
	infile1 = open(infilename1,"r")  # Open secondfile for reading0,
	
	# Read file line by line and write it into a dictionary
	dic={}
	for line0 in infile0:
		paire0 = line0.split()
		dic[paire0[0]] = eval(paire0[1])
		
	for line1 in infile1:
		paire1 = line1.split()
		dic[paire1[0]] = eval(paire1[1])
		
	infile0.close(); infile1.close()
	return dic

if __name__ == '__main__':
    file1 = sys.argv[1]
    file2 = sys.argv[2]	
    mydict = read_file(file1,file2)
    print mydict  


