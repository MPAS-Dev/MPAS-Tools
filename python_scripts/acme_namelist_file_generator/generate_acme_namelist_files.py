#!/usr/bin/env python
import collections
from optparse import OptionParser
import xml.etree.ElementTree as ET

def write_definition_file(registry):#{{{
	definitions = open('namelist_definitions.xml', 'w+')

	# Write definitions header#{{{
	definitions.write('<?xml version="1.0"?>\n')
	definitions.write('\n')
	definitions.write('<namelist_definition>\n')
	definitions.write('\n')
	definitions.write('<!-- Each namelist variable is defined in an <entry> element.  The\n')
	definitions.write('     content of the element is the documentation of how the variable is\n')
	definitions.write('     used.  Other aspects of the variable\'s definition are expressed as\n')
	definitions.write('     attributes of the <entry> element.  Note that it is an XML requirement\n')
	definitions.write('     that the attribute values are enclosed in quotes.  The attributes are:\n')
	definitions.write('\n')
	definitions.write('     id\n')
	definitions.write('          The variable\'s name.  *** N.B. *** The name must be lower case.\n')
	definitions.write('          The module convert all namelist variable names to lower case\n')
	definitions.write('          since Fortran is case insensitive.\n')
	definitions.write('\n')
	definitions.write('     type\n')
	definitions.write('          An abbreviation of the fortran declaration for the variable.\n')
	definitions.write('      Valid declarations are:\n')
	definitions.write('\n')
	definitions.write('          char*n  \n')
	definitions.write('          integer\n')
	definitions.write('          logical\n')
	definitions.write('          real\n')
	definitions.write('\n')
	definitions.write('      Any of these types may be followed by a comma separated list of\n')
	definitions.write('      integers enclosed in parenthesis to indicate an array.\n')
	definitions.write('\n')
	definitions.write('      The current namelist validation code only distinquishes between\n')
	definitions.write('      string and non-string types.\n')
	definitions.write('\n')
	definitions.write('     input_pathname\n')
	definitions.write('          Only include this attribute to indicate that the variable\n')
	definitions.write('          contains the pathname of an input dataset that resides in the\n')
	definitions.write('          CESM inputdata directory tree.  \n')
	definitions.write('\n')
	definitions.write('      The recognized values are "abs" to indicate that an absolute\n')
	definitions.write('          pathname is required, or "rel:var_name" to indicate that the\n')
	definitions.write('          pathname is relative and that the namelist variable "var_name"\n')
	definitions.write('          contains the absolute root directory.\n')
	definitions.write('\n')
	definitions.write('     category\n')
	definitions.write('          A category assigned for organizing the documentation.\n')
	definitions.write('\n')
	definitions.write('     group\n')
	definitions.write('          The namelist group that the variable is declared in.\n')
	definitions.write('\n')
	definitions.write('     valid_values\n')
	definitions.write('          This is an optional attribute that is mainly useful for variables\n')
	definitions.write('          that have only a small number of allowed values.\n')
	definitions.write('                                                                        -->\n')
	definitions.write('\n')
#}}}

	for record in registry.iter("nml_record"):
		record_name = record.attrib["name"]
		definitions.write("\n<!-- %s -->\n"%(record_name.strip()))
		definitions.write("\n")
		for option in record.iter("nml_option"):
			option_name = option.attrib["name"]
			option_type = option.attrib["type"]
			option_description = option.attrib["description"]
			option_units = option.attrib["units"]
			option_default_value = option.attrib["default_value"]
			option_possible_values = option.attrib["possible_values"]

			if ( option_type == "character" ):
				entry_type = "char*1024"
			else:
				entry_type = option_type

			definitions.write("<entry id=\"%s\" type=\"%s\"\n"%(option_name.strip(), entry_type.strip()))
			definitions.write("\tcategory=\"%s\" group=\"%s\">\n"%(record_name.strip(), record_name.strip()))
			definitions.write("%s\n"%(option_description))
			definitions.write("\n")
			definitions.write("Valid values: %s\n"%(option_possible_values))
			definitions.write("Default: %s\n"%(option_default_value))
			definitions.write("</entry>\n\n")


	# Write definitions footer
	definitions.write('</namelist_definition>\n')
	definitions.close()

#}}}

def write_defaults_file(registry, defaults_tree, use_defaults):#{{{
	defaults = open('namelist_defaults.xml', 'w+')

	# Write defaults header#{{{
	defaults.write('<?xml version="1.0"?>\n')
	defaults.write('\n')
	defaults.write('<namelist_defaults>\n')
	defaults.write('\n')
	#}}}

	for record in registry.iter("nml_record"):
		record_name = record.attrib["name"]
		defaults.write("\n<!-- %s -->\n"%(record_name.strip()))
		for option in record.iter("nml_option"):
			option_name = option.attrib["name"]
			option_type = option.attrib["type"]
			option_default_value = option.attrib["default_value"]

			# Extract defaults other than the built in default...
			if use_defaults:
				wrote_opt = False
				for nml_defaults in defaults_tree.iter("namelist_defaults"):
					for def_option in nml_defaults.iter("%s"%(option_name.strip())):
 						if ( not len(def_option.attrib) == 0 ):
							defaults.write("<%s"%(option_name))
							for key, val in def_option.attrib.items():
								defaults.write(" %s=\"%s\""%(key.strip(), val.strip()))
							defaults.write(">%s</%s>\n"%(def_option.text, option_name.strip()))
							wrote_opt = True
						elif len(def_option.attrib) == 0:
							if not option_default_value.strip() == def_option.text.strip().strip("'"):
								print " Default values don't match for option: %s"%(option_name.strip())
								print " Values are: Reg (%s) Def (%s)"%(option_default_value.strip(), def_option.text.strip())
								print " Writing both. Clean up manually..."

								if option_type == 'character':
									defaults.write("<%s>'%s'</%s> <!-- From registry -->\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))
								else:
									defaults.write("<%s>%s</%s> <!-- From registry -->\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))

								defaults.write("<%s>%s</%s> <!-- From old defaults -->\n"%(option_name.strip(), def_option.text.strip(), option_name.strip()))
							else:
								defaults.write("<%s>%s</%s>\n"%(option_name.strip(), def_option.text.strip(), option_name.strip()))
							wrote_opt = True
				if not wrote_opt:
					if option_type == 'character':
						defaults.write("<%s>'%s'</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))
					else:
						defaults.write("<%s>%s</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))
					wrote_opt = True
			else:
				# Write the built in default.
				if option_type == 'character':
					defaults.write("<%s>'%s'</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))
				else:
					defaults.write("<%s>%s</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))


	# Write definitions footer
	defaults.write('\n')
	defaults.write('</namelist_defaults>\n')
	defaults.close()

#}}}

def write_build_namelist_section(registry):#{{{
	build_nml = open('build-namelist-section', 'w+')

	for record in registry.iter("nml_record"):
		record_name = record.attrib["name"]
		message="# Namelist group: %s #\n"%(record_name.strip())

		comment_message=""
		first = True
		for i in message:
			if first:
				comment_message=""
				first = False
			else:
				comment_message = comment_message + "#"
		comment_message = comment_message + "\n"

		build_nml.write(comment_message)
		build_nml.write(message)
		build_nml.write(comment_message)
		build_nml.write("\n")

		for option in record.iter("nml_option"):
			option_name = option.attrib["name"]
			build_nml.write("add_default($nl, '%s');\n"%(option_name.strip()))
	
		build_nml.write("\n")

	build_nml.close()
#}}}

parser = OptionParser()
parser.add_option("-r", "--registry", dest="registry", help="Path to Preprocessed Registry file", metavar="FILE")
parser.add_option("-d", "--defaults", dest="defaults", help="Path to namelist defaults file (Optional)", metavar="FILE")

options, args = parser.parse_args()

if not options.registry:
	parser.error("A processed registry file is required.")

if not options.defaults:
	use_defaults = False
else:
	use_defaults = True

try:
	registry_tree = ET.parse(options.registry)
except:
	parser.error('%s does not exist or is not parsable. Exiting.\nSometimes blank lines at the beginning of the file break the parser.'%options.registry)

if use_defaults:
	try:
		defaults_tree = ET.parse(options.defaults)
	except:
		parser.error('%s does not exist or is not parsable. Exiting.\nSometimes blank lines at the beginning of the file break the parser.'%options.defaults)

	defaults = defaults_tree.getroot()
else:
	defaults = ''

registry = registry_tree.getroot()

write_definition_file(registry)
write_defaults_file(registry, defaults, use_defaults)
write_build_namelist_section(registry)

print "NOTE: You should ensure the build-namelist generated has the correct syntax. Some options need additional logic this script is not capable of generating."
