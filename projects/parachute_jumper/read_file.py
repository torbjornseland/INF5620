#!/usr/bin/env python
import sys

""" This program read data from file and return a dictionary"""

def read_file(infilename0, infilename1):
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
		dic[paire0[0]] = paire0[1]
		
	for line1 in infile1:
		paire1 = line1.split()
		dic[paire1[0]] = paire1[1]
		
	infile0.close(); infile1.close()
	return dic

file1 = sys.argv[1]
file2 = sys.argv[2]	
mydict = read_file(file1,file2)
print mydict  


