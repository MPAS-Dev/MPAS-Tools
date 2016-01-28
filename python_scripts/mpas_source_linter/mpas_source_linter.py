#!/usr/bin/env python
"""
Name: mpas_source_linter.py
Author: Doug Jacobsen
Date: 01/21/2016 

This script can be used to detect syntax errors in source code that may not
cause compilation errors. In general it should help developer enforce a
standard practice of code format while developing their code.

The ouptut from this script is <input_file>.e which is a file containing any
errors detected when parsing the specific file.
"""

import sys, os, glob, shutil, numpy, math
import fnmatch
import argparse
import re
import xml.etree.ElementTree as ET

try:
	from collections import defaultdict
except ImportError:
	from utils import defaultdict

# Define and process input arguments
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--file", dest="file", help="File that should be linted.", metavar="FILE")
args = parser.parse_args()

if not os.path.exists(args.file):
	parser.error('ERROR: File %s does not exist.'%(args.file))

lint_file = open(args.file, 'r')

xml_file = False
fortran_file = False
c_file = False
error_file_open = False
num_errors = 0

# Check for XML Files
if fnmatch.fnmatch(args.file, '*.xml'):
	xml_file = True
	bad_indent = '   '

# Check for Fortran Files
if fnmatch.fnmatch(args.file, '*.F') or fnmatch.fnmatch(args.file, '*.F90'):
	fortran_file = True
	bad_indent = '\t'

if fnmatch.fnmatch(args.file, '*.inc'):
	fortran_file = True
	bad_indent = '\t'

# Check for C/C++ Files
if fnmatch.fnmatch(args.file, '*.c') or fnmatch.fnmatch(args.file, '*.cpp'):
	c_file = True
	bad_indent = '   '

line_num = 0
for block in iter(lambda: lint_file.readline(), ""):
	line_num = line_num + 1
	line_length = len(block)-1
	white_error = False
	len_error = False
	indent_error = False
	kind_error = False

	if ( line_length > 0 ):
		# Test for trailing whitespaces
		match = re.search('[ \t]+$', block)
		if ( not match == None ):
			white_error = True

		# Test for indentation issues
		if ( xml_file or c_file ):
			match = re.search('%s+'%bad_indent, block)
			if ( not match == None ):
				indent_error = True
		elif ( fortran_file ):
			match = re.search('%s+'%bad_indent, block)
			if ( not match == None ):
				indent_error = True

		# Test for line length issues
		if ( fortran_file ) :
			if ( line_length >= 132 ):
				len_error = True

		# Test for issues adding _RKIND to reals
		# Fix to ignore comments, and ignore reals in strings
 		if ( fortran_file ):
			part_arr = block.split('!')
			part = part_arr[0]
			in_str = []

			if ( len(part) > 0 ):
				if part[0] == '"' or part[0] == "'":
					in_str.append(True)
				else:
					in_str.append(False)

				for i in numpy.arange(1, len(part)):
					if part[i] == '"' or part[i] == "'":
						in_str.append(not in_str[-1])
					else:
						in_str.append(in_str[-1])

				starts = []
				ends = []
				
				# Build list of all reals to split on
				for match in re.finditer('[0-9]+\.[0-9]+',part):
					if not in_str[match.start()]:
						starts.append(match.start())
						ends.append(match.end())

				# Test all but the last match
				if ( len(ends) > 1 ):
					for i in numpy.arange(1, len(ends)-2):
						start = ends[i]
						end = starts[i+1]

						match = re.search('_', part[start:end])
						if ( not match ):
							kind_error = True

				# Test last match:
				if ( len(ends) >= 1 ):
					start = ends[-1]
					match = re.search('_', part[start:])
					if ( not match ):
						kind_error = True

				del starts
				del ends
			del in_str

		if ( white_error or len_error or indent_error or kind_error ):

			if not error_file_open:
				error_file = open("%s.e"%(args.file), 'w')
				error_file_open = True

			error_file.write("*************************************************************\n")
			error_file.write("Line number: %d\n"%(line_num))
			error_file.write("Errors:\n")
			if ( indent_error ):
				error_file.write("\tIncorrect characters for indentation\n")
			if ( white_error ):
				error_file.write("\tTrailing white spaces or tabs\n")
			if ( len_error ):
				error_file.write("\tLine is too long (must be less than 132 chars)\n")
				error_file.write("\tLength: %d\n"%(line_length))
			if ( kind_error ):
				error_file.write("\tReal is missing a kind specifier (like _RKIND)\n")
			error_file.write("*************************************************************\n")
			error_file.write("\n")
			num_errors = num_errors + 1

# Attempt to parse an XML file
if ( xml_file ):
 	try:
		tree = ET.parse(args.file)
		root = tree.getroot()
		del root
		del tree
 	except ET.ParseError as e:
		if not error_file_open:
			error_file = open("%s.e"%(args.file), 'w')
			error_file_open = True
		error_file.write("\n")
		error_file.write(" *** ERROR PARSING XML FILE *** \n")
		error_file.write("Error is: %s\n"%(e))
		error_file.write("\n")
		error_file.write(" If this is a full file, the issue may be valid. If it's part of a file, you can probably ignore this error.\n")

lint_file.close()
if error_file_open:
	error_file.close()
sys.exit(num_errors)
