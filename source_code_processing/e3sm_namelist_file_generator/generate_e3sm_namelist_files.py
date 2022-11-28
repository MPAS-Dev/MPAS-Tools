#!/usr/bin/env python
"""
This script can be used to autogenerate some of the files used by ACME/CESM to
define namelist options for each component.  These files are:

* models/<COMPONENT>/bld/build-namelist:
  This perl script builds the namelist that ACME uses for the component.
  It creates the default namelist that the component will use by a hardcoded
  list of namelist options that should be included.
  It uses the namelist_definition_*.xml file to validate the data type of
  each namelist option that it attempts to add to the default namelist.  It then
  uses namelist_defaults_*.xml  to assign default values for each of these options.

* models/<COMPONENT>/bld/namelist-files/namelist_definition_<MODEL>.xml:
  This xml file defines the namelist options that exist for this component,
  including a description of each that is used for documentation purposes.

* models/<COMPONENT>/bld/namelist-files/namelist_defaults_<MODEL>.xml:
  This xml file defines the default value for each namelist option defined in
  namelist_definition_*.xml.  The default values defined here for ACME could differ
  from the defaults used in standalone MPAS.

This script requires a processed registry file from a compiled instance of MPAS
(e.g., src/core_<CORE>/Registry_processed.xml).  The Registry definitions of MPAS
namelist options are used to automatically generate:
* namelist_definitions.xml: the entire namelist_definition_*.xml file
* namelist_defaults.xml: the entire namelist_defaults_*.xml file
* build-namelist-group-list: the section of the build-namelist file that defines
  which namelist groups exist.  The content of the file created should be manually
  copied into build-namelist in the appropriate location, replacing its previous contents.
* build-namelist-section:  the section of the build-namelist file that adds
  namelist options to the default  namelist file.  Every namelist option found in
  the Registry file  specified will be added.  The content of the file created
  should then be manually copied into  build-namelist in the appropriate location,
  replacing its previous contents.

Optionally, a pre-existing namelist_defaults_*.xml file can be specified.  If it is,
there will be warnings to the screen and lines written to the newly generated
namelist_default_*.xml file that is output indicating where the values in the
pre-existing namelist_defaults file differ from those in the MPAS Registry file.
This is useful when the defaults used for ACME differ from those used in standalone
MPAS, and the ACME defaults should be preserved.  This information alerts the operator
to these differences and allows the conflicts to be manually resolved.

"""
from __future__ import print_function
import argparse
import xml.etree.ElementTree as ET


def write_definition_file(registry):#{{{
	nl_defin_filename = "namelist_definitions.xml"
	print("=== Writing namelist definitions to: " + nl_defin_filename)
	definitions = open(nl_defin_filename, 'w+')

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
	definitions.write('          char*n\n')
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
	definitions.write('          CESM inputdata directory tree.\n')
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
			try:
				option_description = option.attrib["description"]
			except:
				option_description = 'MISSING DESCRIPTION'

			try:
				option_units = option.attrib["units"]
			except:
				option_units = 'MISSING UNITS'

			option_default_value = option.attrib["default_value"]

			try:
				option_possible_values = option.attrib["possible_values"]
			except:
				option_possible_values = 'MISSING POSSIBLE VALUES'

			if ( option_type == "character" ):
				entry_type = "char*1024"
			else:
				entry_type = option_type

			definitions.write("<entry id=\"%s\" type=\"%s\"\n"%(option_name.strip(), entry_type.strip()))
			definitions.write("\tcategory=\"%s\" group=\"%s\">\n"%(record_name.strip(), record_name.strip()))
			definitions.write("%s\n"%(option_description))
			definitions.write("\n")
			if (len(option_possible_values) > 0):
				definitions.write("Valid values: %s\n"%(option_possible_values))
			else:
				definitions.write("Valid values:\n")
			definitions.write("Default: Defined in namelist_defaults.xml\n")
			definitions.write("</entry>\n\n")


	# Write definitions footer
	definitions.write('</namelist_definition>\n')
	definitions.close()
	print("=== Complete.\n")

#}}}

def write_defaults_file(registry, defaults_tree, use_defaults):#{{{
	nl_defaults_filename = "namelist_defaults.xml"
	print("=== Writing namelist defaults to: " + nl_defaults_filename)
	defaults = open(nl_defaults_filename, 'w+')

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
								print("  Default values don't match for option: %s"%(option_name.strip()))
								print("    Values are: Reg (%s) Def (%s)"%(option_default_value.strip(), def_option.text.strip()))
								print("    Writing both. Clean up manually...")

								if option_type == 'character':
									defaults.write("<<<<<<<<< FROM REGISTRY\n")
									defaults.write("<%s>'%s'</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))
								else:
									defaults.write("<<<<<<<<< FROM REGISTRY\n")
									defaults.write("<%s>%s</%s>\n"%(option_name.strip(), option_default_value.strip(), option_name.strip()))

								defaults.write("=======================\n")
								defaults.write("<%s>%s</%s>\n"%(option_name.strip(), def_option.text.strip(), option_name.strip()))
								defaults.write("<<<<<<<<< FROM OLD DEFAULTS\n")
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
	print("=== Complete.\n")

#}}}

def write_build_namelist_section(registry):#{{{
	bld_nl_sec_filename = "build-namelist-section-new"
	print("=== Writing section of the 'build-namelist' file that writes a list of namelist options to: " + bld_nl_sec_filename)
	build_nml = open(bld_nl_sec_filename, 'w+')

	group_list = open('build-namelist-group-list', 'w+')

	group_list.write('my @groups = qw(')

	for record in registry.iter("nml_record"):
		record_name = record.attrib["name"]
		message="# Namelist group: %s #\n"%(record_name.strip())
		group_list.write('%s\n                '%(record_name.strip().lower()))

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

	group_list.write(');\n')
	group_list.close()
	build_nml.close()
	print("   NOTE: You should ensure the build-namelist generated has the correct syntax. Some options need additional logic this script is not capable of generating.")
	print("=== Complete.\n")
#}}}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-r", "--registry", dest="registry", help="Path to Preprocessed Registry file", metavar="FILE", required=True)
parser.add_argument("-d", "--defaults", dest="defaults", help="Path to namelist defaults file (Optional)", metavar="FILE")

args = parser.parse_args()

if not args.registry:
	parser.error("A processed registry file is required.")

if not args.defaults:
	use_defaults = False
else:
	use_defaults = True

registry_tree = ET.parse(args.registry)

if use_defaults:
	defaults_tree = ET.parse(args.defaults)

	defaults = defaults_tree.getroot()
else:
	defaults = ''

registry = registry_tree.getroot()

write_definition_file(registry)
write_defaults_file(registry, defaults, use_defaults)
write_build_namelist_section(registry)

