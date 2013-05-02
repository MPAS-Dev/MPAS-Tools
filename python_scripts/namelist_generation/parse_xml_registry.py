#!/usr/bin/python
from collections import defaultdict
from optparse import OptionParser
import xml.etree.ElementTree as ET

parser = OptionParser()
parser.add_option("-f", "--file", dest="registry_path", help="Path to Registry file", metavar="FILE")
parser.add_option("-d", "--tex_dir", dest="latex_dir", help="Path to directory with latex addition files.", metavar="DIR")
parser.add_option("-p", "--tex_path", dest="latex_path", help="Path to latex input files that will be written to generated latex.", metavar="PATH")

options, args = parser.parse_args()

if not options.registry_path:
	parser.error("Registry file is required")

if not options.latex_dir:
	print 'Directory with group latex files is missing. Skipping addition of latex files.'
	extra_latex = False
else:
	if not options.latex_path:
		parser.error('Need latex path with latex directory.')
	extra_latex = True


latex_missing_string = '{\\bf \color{red} MISSING}'
dimension_table_header = '{\\bf Name} & {\\bf Units} & {\\bf Description}'
variable_table_header = '{\\bf Name} & {\\bf Description}'
namelist_table_header = '{\\bf Name} & {\\bf Description}'

registry_path = options.registry_path

try:
	registry_tree = ET.parse(registry_path)
except:
	parser.error('%s does not exist or is not parsable. Exiting.'%registry_path)

registry = registry_tree.getroot()

# Write default namelist
namelist = open('namelist.input.generated', 'w+')
for nml_rec in registry.iter("nml_record"):
	namelist.write('&%s\n'%nml_rec.attrib['name'])
	for nml_opt in nml_rec.iter("nml_option"):
		if nml_opt.attrib['type'] == "character":
			namelist.write('\t%s = "%s"\n'%(nml_opt.attrib['name'], nml_opt.attrib['default_value']))
		else:
			namelist.write('\t%s = %s\n'%(nml_opt.attrib['name'], nml_opt.attrib['default_value']))

	namelist.write('/\n')

# Write file that defines version string for model.
latex = open('define_version.tex', 'w+')
try:
	version_string = registry.attrib['version']
except:
	version_string = '{\\bf MISSING}'
latex.write('\\newcommand{\\version}{%s}\n'%version_string)
latex.close()

# Write dimension table documentation latex file.
latex = open('dimension_table_documentation.tex', 'w+')
latex.write('\chapter{Dimensions}\n')
latex.write('\label{chap:dimensions}\n')
latex.write('{\small\n')
latex.write('\\begin{center}\n')
latex.write('\\begin{longtable}{| p{1.0in} || p{1.0in} | p{4.0in} |}\n')
latex.write('	\hline \n')
latex.write('	%s \\\\\n'%dimension_table_header)
latex.write('	\hline \n')
latex.write('	\hline \n')
for dims in registry.iter("dims"):
	for dim in dims.iter("dim"):
		dim_name = dim.attrib['name']
		try:
			dim_description = dim.attrib['description']
		except:
			dim_description = latex_missing_string

		try:
			dim_units = dim.attrib['units']
			if dim_units == "":
				dim_units = latex_missing_string
			else:
				dim_units = "$%s$"%dim_units.replace(' ', '$ $')
		except:
			dim_units = latex_missing_string

		if dim_description == "":
			dim_description = latex_missing_string
		else:
			equations = dim_description.find('$')
			if equations != -1:
				dim_desc_split = dim_description.split("$")

				if dim_description.replace('_','')[0] == "$":
					replace = False
					dim_description = "$"
				else:
					replace = True
					dim_description = ""

				for part in dim_desc_split:
					if replace:
						dim_description = "%s %s"%(dim_description, part.replace('_','\_'))
						replace = False
					else:
						dim_description = "%s $%s$"%(dim_description, part)
						replace = True
			else:
				dim_description = "%s"%dim_description.replace('_','\_')

		latex.write('	%s & %s & %s \\\\ \n'%(dim_name.replace('_','\_'), dim_units.replace('_','\_'), dim_description.replace('_','\_')))
		latex.write('	\hline\n')

latex.write('\end{longtable}\n')
latex.write('\end{center}\n')
latex.write('}\n')
latex.close()


# Write namelist table documentation latex file.
latex = open('namelist_table_documentation.tex', 'w+')
latex.write('\chapter[Namelist options]{\hyperref[chap:namelist_sections]{Namelist options}}\n')
latex.write('\label{chap:namelist_tables}\n')
latex.write('Embedded links point to more detailed namelist information in the appendix.\n')
for nml_rec in registry.iter("nml_record"):
	rec_name = nml_rec.attrib['name']
	latex.write('\section[%s]{\hyperref[sec:nm_sec_%s]{%s}}\n'%(rec_name.replace('_','\_'), rec_name, rec_name.replace('_','\_')))
	latex.write('\label{sec:nm_tab_%s}\n'%(rec_name))

	# Add input line if file exists.
	try:
		junk_file = open('%s/%s.tex'%(options.latex_dir,rec_name), 'r')
		latex.write('\input{%s/%s.tex}\n'%(options.latex_path, rec_name))
		junk_file.close()
	except:
		latex.write('')

	latex.write('{\small\n')
	latex.write('\\begin{center}\n')
	latex.write('\\begin{longtable}{| p{2.0in} || p{4.0in} |}\n')
	latex.write('	\hline\n')
	latex.write('	%s \\\\\n'%namelist_table_header)
	latex.write('	\hline\n')
	latex.write('	\hline\n')

	for nml_opt in nml_rec.iter("nml_option"):
		opt_name = nml_opt.attrib['name']

		try:
			opt_description = nml_opt.attrib['description']
		except:
			opt_description = latex_missing_string

		if opt_description == "":
			opt_description = latex_missing_string
		else:
			equations = opt_description.find('$')
			if equations != -1:
				opt_desc_split = opt_description.split("$")

				if opt_description.replace(' ','')[0] == "$":
					replace = False
					opt_description = "$"
				else:
					replace = True
					opt_description = ""

				for part in opt_desc_split:
					if replace:
						opt_description = "%s %s"%(opt_description, part.replace('_','\_'))
						replace = False
					else:
						opt_description = "%s $%s$"%(opt_description, part)
						replace = True
			else:
				opt_description = "%s"%opt_description.replace('_','\_')

		latex.write('	\hyperref[subsec:nm_sec_%s]{%s} & %s \\\\\n'%(opt_name, opt_name.replace('_','\_'), opt_description))
		latex.write('	\hline\n')

	latex.write('\end{longtable}\n')
	latex.write('\end{center}\n')
	latex.write('}\n')
latex.close()
		
# Write namelist section documentation latex file.
latex = open('namelist_section_documentation.tex', 'w+')
latex.write('\chapter[Namelist options]{\hyperref[chap:namelist_tables]{Namelist options}}\n')
latex.write('\label{chap:namelist_sections}\n')
latex.write('Embedded links point to information in chapter \\ref{chap:namelist_tables}\n')
for nml_rec in registry.iter("nml_record"):
	rec_name = nml_rec.attrib["name"]

	latex.write('\section[%s]{\hyperref[sec:nm_tab_%s]{%s}}\n'%(rec_name.replace('_','\_'),rec_name,rec_name.replace('_','\_')))
	latex.write('\label{sec:nm_sec_%s}\n'%rec_name)

	for nml_opt in nml_rec.iter("nml_option"):
		opt_name = nml_opt.attrib["name"]
		opt_type = nml_opt.attrib["type"]
		opt_value = nml_opt.attrib["default_value"]

		try:
			opt_possible_values = nml_opt.attrib["possible_values"]
		except:
			opt_possible_values = latex_missing_string

		try:
			opt_units = nml_opt.attrib["units"]
			if opt_units == "":
				opt_units = latex_missing_string
			else:
				opt_units = "$%s$"%opt_units.replace(' ', '$ $')
		except:
			opt_units = latex_missing_string

		try:
			opt_description = nml_opt.attrib["description"]
		except:
			opt_description = latex_missing_string

		if opt_possible_values == "":
			opt_possible_values = latex_missing_string


		if opt_description == "":
			opt_description = latex_missing_string.replace('_','\_')
		else:
			equations = opt_description.find('$')
			if equations != -1:
				opt_desc_split = opt_description.split("$")

				if opt_description.replace('_','')[0] == "$":
					replace = False
					opt_description = "$"
				else:
					replace = True
					opt_description = ""

				for part in opt_desc_split:
					if replace:
						opt_description = "%s %s"%(opt_description, part.replace('_','\_'))
						replace = False
					else:
						opt_description = "%s $%s$"%(opt_description, part)
						replace = True
			else:
				opt_description = "%s"%opt_description.replace('_','\_')

		latex.write('\subsection[%s]{\hyperref[sec:nm_tab_%s]{%s}}\n'%(opt_name.replace('_','\_'),rec_name,opt_name.replace('_','\_')))
		latex.write('\label{subsec:nm_sec_%s}\n'%opt_name)
		latex.write('\\begin{center}\n')
		latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
		latex.write('    \hline\n')
		latex.write('    Type: & %s \\\\\n'%opt_type.replace('_','\_'))
		latex.write('    \hline\n')
		latex.write('    Units: & %s \\\\\n'%opt_units.replace('_','\_'))
		latex.write('    \hline\n')
		latex.write('    Default Value: & %s \\\\\n'%opt_value.replace('_','\_'))
		latex.write('    \hline\n')
		latex.write('    Possible Values: & %s \\\\\n'%opt_possible_values.replace('_','\_'))
		latex.write('    \hline\n')
		latex.write('    \caption{%s: %s}\n'%(opt_name.replace('_','\_'), opt_description))
		latex.write('\end{longtable}\n')
		latex.write('\end{center}\n')
latex.close()

# Write variable table documentation latex file
latex = open('variable_table_documentation.tex', 'w+')
latex.write('\chapter[Variable definitions]{\hyperref[chap:variable_sections]{Variable definitions}}\n')
latex.write('\label{chap:variable_tables}\n')
latex.write('Embedded links point to more detailed variable information in the appendix.\n')
for var_struct in registry.iter("var_struct"):
	struct_name = var_struct.attrib['name']
	latex.write('\section[%s]{\hyperref[sec:var_sec_%s]{%s}}\n'%(struct_name.replace('_','\_'),struct_name,struct_name.replace('_','\_')))
	latex.write('\label{sec:var_tab_%s}\n'%struct_name)

	try:
		junk_file = open('%s/%s.tex'%(options.latex_dir,struct_name), 'r')
		latex.write('\input{%s/%s.tex}\n'%(options.latex_path, struct_name))
		junk_file.close()
	except:
		latex.write('')

	latex.write('{\small\n')
	latex.write('\\begin{center}\n')
	latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
	latex.write('	\hline\n')
	latex.write('	%s \\\\\n'%variable_table_header)
	latex.write('	\hline\n')

	for node in var_struct.getchildren():
		if node.tag == 'var_array':
			for var in node.iter("var"):
				var_name = var.attrib['name']
				var_description = var.attrib['description']

				if var_description == "":
					var_description = latex_missing_string.replace('_','\_')
				else:
					equations = var_description.find('$')
					if equations != -1:
						var_desc_split = var_description.split("$")

						if var_description.replace('_','')[0] == "$":
							replace = False
							var_description = "$"
						else:
							replace = True
							var_description = ""

						for part in var_desc_split:
							if replace:
								var_description = "%s %s"%(var_description, part.replace('_','\_'))
								replace = False
							else:
								var_description = "%s $%s$"%(var_description, part)
								replace = True
					else:
						var_description = "%s"%var_description.replace('_','\_')

				latex.write('	\hyperref[subsec:var_sec_%s_%s]{%s} & %s \\\\\n'%(struct_name, var_name, var_name.replace('_','\_'), var_description))
				latex.write('	\hline\n')
		elif node.tag == 'var':
			var = node
			var_name = var.attrib['name']
			var_description = var.attrib['description']

			if var_description == "":
				var_description = latex_missing_string.replace('_','\_')
			else:
				equations = var_description.find('$')
				if equations != -1:
					var_desc_split = var_description.split("$")

					if var_description.replace('_','')[0] == "$":
						replace = False
						var_description = "$"
					else:
						replace = True
						var_description = ""

					for part in var_desc_split:
						if replace:
							var_description = "%s %s"%(var_description, part.replace('_','\_'))
							replace = False
						else:
							var_description = "%s $%s$"%(var_description, part)
							replace = True
				else:
					var_description = "%s"%var_description.replace('_','\_')

			latex.write('	\hyperref[subsec:var_sec_%s_%s]{%s} & %s \\\\\n'%(struct_name, var_name, var_name.replace('_','\_'), var_description))
			latex.write('	\hline\n')

	latex.write('\end{longtable}\n')
	latex.write('\end{center}\n')
	latex.write('}\n')
latex.close()

# Write variable section documentation latex file
latex = open('variable_section_documentation.tex', 'w+')
latex.write('\chapter[Variable definitions]{\hyperref[chap:variable_tables]{Variable definitions}}\n')
latex.write('\label{chap:variable_sections}\n')
latex.write('Embedded links point to information in chapter \\ref{chap:variable_tables}\n')
for var_struct in registry.iter("var_struct"):
	struct_name = var_struct.attrib['name']
	struct_time_levs = var_struct.attrib['time_levs']
	latex.write('\section[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(struct_name.replace('_','\_'),struct_name, struct_name.replace('_','\_')))
	latex.write('\label{sec:var_sec_%s}\n'%struct_name)
	
	for node in var_struct.getchildren():
		if node.tag == 'var_array':
			var_arr = node
			for var in var_arr.iter("var"):
				var_arr_name = var_arr.attrib['name']
				var_arr_type = var_arr.attrib['type']
				var_arr_dims = var_arr.attrib['dimensions']

				var_name = var.attrib['name']
				var_arr_group = var.attrib['array_group']

				try:
					var_persistence = var_arr.attrib['persistence']
				except:
					var_persistence = 'persistent'

				try:
					var_name_in_code = var.attrib['name_in_code']
				except:
					var_name_in_code = var_name

				try:
					var_units = var.attrib['units']
					if var_units == "":
						var_units = latex_missing_string
					else:
						var_units = "$%s$"%var_units.replace(' ', '$ $')
				except:
					var_units = latex_missing_string

				try:
					var_description = var.attrib['description']
				except:
					var_description = latex_missing_string.replace('_','\_')

				try:
					var_streams = var.attrib['streams'].replace('s','Sfc ').replace('i','Input ').replace('r', 'Restart ').replace('o','Output ')
				except:
					var_streams = "None"

				if int(struct_time_levs) > 1:
					var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code)
					var_path = "domain %% blocklist %% %s %% time_levs(:) %% %s %% %s"%(struct_name, struct_name, var_arr_name)
				else:
					var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code)
					var_path = "domain %% blocklist %% %s %% %s"%(struct_name, var_arr_name)

				if var_description == "":
					var_description = latex_missing_string.replace('_','\_')
				else:
					equations = var_description.find('$')
					if equations != -1:
						var_desc_split = var_description.split("$")

						if var_description.replace('_','')[0] == "$":
							replace = False
							var_description = "$"
						else:
							replace = True
							var_description = ""

						for part in var_desc_split:
							if replace:
								var_description = "%s %s"%(var_description, part.replace('_','\_'))
								replace = False
							else:
								var_description = "%s $%s$"%(var_description, part)
								replace = True
					else:
						var_description = "%s"%var_description.replace('_','\_')

				
				latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(var_name.replace('_','\_'),struct_name, var_name.replace('_','\_')))
				latex.write('\label{subsec:var_sec_%s_%s}\n'%(struct_name,var_name))
				# Tabular Format:
				latex.write('\\begin{center}\n')
				latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
				latex.write('        \hline \n')
				latex.write('        Type: & %s \\\\\n'%var_arr_type)
				latex.write('        \hline \n')
				latex.write('        Units: & %s \\\\\n'%var_units)
				latex.write('        \hline \n')
				latex.write('        Dimension: & %s \\\\\n'%var_arr_dims)
				latex.write('        \hline \n')
				latex.write('        Persistence: & %s \\\\\n'%var_persistence)
				latex.write('        \hline \n')
				latex.write('		 Default Streams: & %s \\\\\n'%var_streams)
				latex.write('        \hline \n')
				latex.write('		 Index in %s Array: & %s \\\\\n'%(var_arr_name.replace('_','\_'), var_index.replace('_','\_').replace('%','\%')))
				latex.write('		 \hline \n')
				latex.write('		 Location in code: & %s \\\\\n'%var_path.replace('_','\_').replace('%','\%'))
				latex.write('		 \hline \n')
				latex.write('		 Array Group: & %s \\\\\n'%var_arr_group)
				latex.write('		 \hline \n')
				latex.write('    \caption{%s: %s}\n'%(var_name.replace('_','\_'),var_description))
				latex.write('\end{longtable}\n')
				latex.write('\end{center}\n')
		elif node.tag == 'var':
			var = node
			try:
				# Skip super array variables. They have different tables.
				var_arr_group = var.attrib['array_group']
			except:
				var_name = var.attrib['name']
				var_type = var.attrib['type']
				var_dims = var.attrib['dimensions']

				try:
					var_persistence = var_arr.attrib['persistence']
				except:
					var_persistence = 'persistent'

				try:
					var_name_in_code = var.attrib['name_in_code']
				except:
					var_name_in_code = var_name

				try:
					var_units = var.attrib['units']
					if var_units == "":
						var_units = latex_missing_string
					else:
						var_units = "$%s$"%var_units.replace(' ', '$ $')
				except:
					var_units = latex_missing_string

				try:
					var_description = var.attrib['description']
				except:
					var_description = latex_missing_string.replace('_','\_')

				try:
					var_streams = var.attrib['streams'].replace('s','Sfc ').replace('i','Input ').replace('r', 'Restart ').replace('o','Output ')
				except:
					var_streams = "None"

				if int(struct_time_levs) > 1:
					var_path = "domain %% blocklist %% %s %% time_levs(:) %% %s %% %s"%(struct_name, struct_name, var_name_in_code)
				else:
					var_path = "domain %% blocklist %% %s %% %s"%(struct_name, var_name_in_code)

				if var_description == "":
					var_description = latex_missing_string.replace('_','\_')
				else:
					equations = var_description.find('$')
					if equations != -1:
						var_desc_split = var_description.split("$")

						if var_description.replace('_','')[0] == "$":
							replace = False
							var_description = "$"
						else:
							replace = True
							var_description = ""

						for part in var_desc_split:
							if replace:
								var_description = "%s %s"%(var_description, part.replace('_','\_'))
								replace = False
							else:
								var_description = "%s $%s$"%(var_description, part)
								replace = True
					else:
						var_description = "%s"%var_description.replace('_','\_')

				
				latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(var_name.replace('_','\_'),struct_name, var_name.replace('_','\_')))
				latex.write('\label{subsec:var_sec_%s_%s}\n'%(struct_name,var_name))
				# Tabular Format:
				latex.write('\\begin{center}\n')
				latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
				latex.write('        \hline \n')
				latex.write('        Type: & %s \\\\\n'%var_type)
				latex.write('        \hline \n')
				latex.write('        Units: & %s \\\\\n'%var_units)
				latex.write('        \hline \n')
				latex.write('        Dimension: & %s \\\\\n'%var_dims)
				latex.write('        \hline \n')
				latex.write('        Persistence: & %s \\\\\n'%var_persistence)
				latex.write('        \hline \n')
				latex.write('		 Default Streams: & %s \\\\\n'%var_streams)
				latex.write('        \hline \n')
				latex.write('		 Location in code: & %s \\\\\n'%var_path.replace('_','\_').replace('%','\%'))
				latex.write('		 \hline \n')
				latex.write('    \caption{%s: %s}\n'%(var_name.replace('_','\_'),var_description))
				latex.write('\end{longtable}\n')
				latex.write('\end{center}\n')

#	# Extract variables from variable arrays
#	for var_arr in var_struct.iter("var_array"):
#		for var in var_arr.iter("var"):
#			var_arr_name = var_arr.attrib['name']
#			var_arr_type = var_arr.attrib['type']
#			var_arr_dims = var_arr.attrib['dimensions']
#
#			var_name = var.attrib['name']
#			var_arr_group = var.attrib['array_group']
#
#			try:
#				var_persistence = var_arr.attrib['persistence']
#			except:
#				var_persistence = 'persistent'
#
#			try:
#				var_name_in_code = var.attrib['name_in_code']
#			except:
#				var_name_in_code = var_name
#
#			try:
#				var_units = var.attrib['units']
#				if var_units == "":
#					var_units = latex_missing_string
#				else:
#					var_units = "$%s$"%var_units.replace(' ', '$ $')
#			except:
#				var_units = latex_missing_string
#
#			try:
#				var_description = var.attrib['description']
#			except:
#				var_description = latex_missing_string.replace('_','\_')
#
#			try:
#				var_streams = var.attrib['streams'].replace('s','Sfc ').replace('i','Input ').replace('r', 'Restart ').replace('o','Output ')
#			except:
#				var_streams = "None"
#
#			if int(struct_time_levs) > 1:
#				var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code)
#				var_path = "domain %% blocklist %% %s %% time_levs(:) %% %s %% %s"%(struct_name, struct_name, var_arr_name)
#			else:
#				var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code)
#				var_path = "domain %% blocklist %% %s %% %s"%(struct_name, var_arr_name)
#
#			if var_description == "":
#				var_description = latex_missing_string.replace('_','\_')
#			else:
#				equations = var_description.find('$')
#				if equations != -1:
#					var_desc_split = var_description.split("$")
#
#					if var_description.replace('_','')[0] == "$":
#						replace = False
#						var_description = "$"
#					else:
#						replace = True
#						var_description = ""
#
#					for part in var_desc_split:
#						if replace:
#							var_description = "%s %s"%(var_description, part.replace('_','\_'))
#							replace = False
#						else:
#							var_description = "%s $%s$"%(var_description, part)
#							replace = True
#				else:
#					var_description = "%s"%var_description.replace('_','\_')
#
#			
#			latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(var_name.replace('_','\_'),struct_name, var_name.replace('_','\_')))
#			latex.write('\label{subsec:var_sec_%s_%s}\n'%(struct_name,var_name))
#			# Tabular Format:
#			latex.write('\\begin{center}\n')
#			latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
#			latex.write('        \hline \n')
#			latex.write('        Type: & %s \\\\\n'%var_arr_type)
#			latex.write('        \hline \n')
#			latex.write('        Units: & %s \\\\\n'%var_units)
#			latex.write('        \hline \n')
#			latex.write('        Dimension: & %s \\\\\n'%var_arr_dims)
#			latex.write('        \hline \n')
#			latex.write('        Persistence: & %s \\\\\n'%var_persistence)
#			latex.write('        \hline \n')
#			latex.write('		 Default Streams: & %s \\\\\n'%var_streams)
#			latex.write('        \hline \n')
#			latex.write('		 Index in %s Array: & %s \\\\\n'%(var_arr_name.replace('_','\_'), var_index.replace('_','\_').replace('%','\%')))
#			latex.write('		 \hline \n')
#			latex.write('		 Location in code: & %s \\\\\n'%var_path.replace('_','\_').replace('%','\%'))
#			latex.write('		 \hline \n')
#			latex.write('		 Array Group: & %s \\\\\n'%var_arr_group)
#			latex.write('		 \hline \n')
#			latex.write('    \caption{%s: %s}\n'%(var_name.replace('_','\_'),var_description))
#			latex.write('\end{longtable}\n')
#			latex.write('\end{center}\n')
#
#	# Parse variables not in variable arrays
#	for var in var_struct.iter("var"):
#		try:
#			# Skip super array variables. They have different tables.
#			var_arr_group = var.attrib['array_group']
#		except:
#			var_name = var.attrib['name']
#			var_type = var.attrib['type']
#			var_dims = var.attrib['dimensions']
#
#			try:
#				var_persistence = var_arr.attrib['persistence']
#			except:
#				var_persistence = 'persistent'
#
#			try:
#				var_name_in_code = var.attrib['name_in_code']
#			except:
#				var_name_in_code = var_name
#
#			try:
#				var_units = var.attrib['units']
#				if var_units == "":
#					var_units = latex_missing_string
#				else:
#					var_units = "$%s$"%var_units.replace(' ', '$ $')
#			except:
#				var_units = latex_missing_string
#
#			try:
#				var_description = var.attrib['description']
#			except:
#				var_description = latex_missing_string.replace('_','\_')
#
#			try:
#				var_streams = var.attrib['streams'].replace('s','Sfc ').replace('i','Input ').replace('r', 'Restart ').replace('o','Output ')
#			except:
#				var_streams = "None"
#
#			if int(struct_time_levs) > 1:
#				var_path = "domain %% blocklist %% %s %% time_levs(:) %% %s %% %s"%(struct_name, struct_name, var_name_in_code)
#			else:
#				var_path = "domain %% blocklist %% %s %% %s"%(struct_name, var_name_in_code)
#
#			if var_description == "":
#				var_description = latex_missing_string.replace('_','\_')
#			else:
#				equations = var_description.find('$')
#				if equations != -1:
#					var_desc_split = var_description.split("$")
#
#					if var_description.replace('_','')[0] == "$":
#						replace = False
#						var_description = "$"
#					else:
#						replace = True
#						var_description = ""
#
#					for part in var_desc_split:
#						if replace:
#							var_description = "%s %s"%(var_description, part.replace('_','\_'))
#							replace = False
#						else:
#							var_description = "%s $%s$"%(var_description, part)
#							replace = True
#				else:
#					var_description = "%s"%var_description.replace('_','\_')
#
#			
#			latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(var_name.replace('_','\_'),struct_name, var_name.replace('_','\_')))
#			latex.write('\label{subsec:var_sec_%s_%s}\n'%(struct_name,var_name))
#			# Tabular Format:
#			latex.write('\\begin{center}\n')
#			latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
#			latex.write('        \hline \n')
#			latex.write('        Type: & %s \\\\\n'%var_type)
#			latex.write('        \hline \n')
#			latex.write('        Units: & %s \\\\\n'%var_units)
#			latex.write('        \hline \n')
#			latex.write('        Dimension: & %s \\\\\n'%var_dims)
#			latex.write('        \hline \n')
#			latex.write('        Persistence: & %s \\\\\n'%var_persistence)
#			latex.write('        \hline \n')
#			latex.write('		 Default Streams: & %s \\\\\n'%var_streams)
#			latex.write('        \hline \n')
#			latex.write('		 Location in code: & %s \\\\\n'%var_path.replace('_','\_').replace('%','\%'))
#			latex.write('		 \hline \n')
#			latex.write('    \caption{%s: %s}\n'%(var_name.replace('_','\_'),var_description))
#			latex.write('\end{longtable}\n')
#			latex.write('\end{center}\n')
latex.close()




