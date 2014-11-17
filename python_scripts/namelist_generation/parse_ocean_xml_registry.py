#!/usr/bin/env python
import collections
from optparse import OptionParser
import xml.etree.ElementTree as ET

parser = OptionParser()
parser.add_option("-f", "--forward_registry", dest="forward_registry_path", help="Path to Mode Forward's Preprocessed Registry file", metavar="FILE")
parser.add_option("-a", "--analysis_registry", dest="analysis_registry_path", help="Path to Mode Analysis' Preprocessed Registry file", metavar="FILE")
parser.add_option("-d", "--tex_dir", dest="latex_dir", help="Path to directory with latex addition files.", metavar="DIR")
parser.add_option("-p", "--tex_path", dest="latex_path", help="Path to latex input files that will be written to generated latex.", metavar="PATH")

options, args = parser.parse_args()

def break_string(string):#{{{
	i = 0.0
	idx = -1

	size = 0

	big_size = 1.8
	small_size = 1.2
	really_small_size = 0.2

	big_count = 0
	small_count = 0
	really_small_count = 0

	for c in string:
		idx = idx + 1
		if c.isupper():
			big_count = big_count + 1
			size = size + big_size
		else:
			if c == "l" or c == "i":
				really_small_count = really_small_count + 1
				size = size + really_small_size
			else:
				small_count = small_count + 1
				size = size + small_size

		if size >= 33.5:
			return idx

	return -1
	#}}}

def write_dimension_table(latex, registry, mode):#{{{
	latex_missing_string = '{\\bf \color{red} MISSING}'
	dimension_table_header = '{\\bf Name} & {\\bf Units} & {\\bf Description}'

	latex.write('\section{Dimensions}\n')
	latex.write('\label{sec:%s_dimensions}\n'%(mode))
	latex.write('{\small\n')
	latex.write('\\begin{center}\n')
	latex.write('\\begin{longtable}{| p{1.2in} || p{1.0in} | p{4.0in} |}\n')
	latex.write('	\hline \n')
	latex.write('	%s \\endfirsthead\n'%dimension_table_header)
	latex.write('	\hline \n')
	latex.write('	%s (Continued) \\endhead\n'%dimension_table_header)
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
#}}}

def write_namelist_table(latex, registry, mode):#{{{
	latex_missing_string = '{\\bf \color{red} MISSING}'
	namelist_table_header = '{\\bf Name} & {\\bf Description}'

	latex.write('\section[Namelist options]{\hyperref[chap:namelist_sections]{Namelist options}}\n')
	latex.write('\label{sec:%s_namelist_tables}\n'%(mode))
	latex.write('Embedded links point to more detailed namelist information in the appendix.\n')
	for nml_rec in registry.iter("nml_record"):
		rec_name = nml_rec.attrib['name']
		#latex.write('\subsection[%s]{\hyperref[sec:nm_sec_%s]{%s}}\n'%(rec_name.replace('_','\_'), rec_name, rec_name.replace('_','\_')))
		latex.write('\subsection[%s]{%s}\n'%(rec_name.replace('_','\_'), rec_name.replace('_','\_')))
		latex.write('\label{subsec:%s_nm_tab_%s}\n'%(mode, rec_name))

		# Add input line if file exists.
		try:
			junk_file = open('%s/%s.tex'%(options.latex_dir,rec_name), 'r')
			latex.write('\input{%s/%s.tex}\n'%(options.latex_path, rec_name))
			junk_file.close()
		except:
			latex.write('')

		latex.write('\\vspace{0.5in}\n')
		latex.write('{\small\n')
		latex.write('\\begin{center}\n')
		latex.write('\\begin{longtable}{| p{2.0in} || p{4.0in} |}\n')
		latex.write('	\hline\n')
		latex.write('	%s \\endfirsthead\n'%namelist_table_header)
		latex.write('	\hline \n')
		latex.write('	%s (Continued) \\endhead\n'%namelist_table_header)
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

			idx = break_string(opt_name)
			if idx >= 29:
				latex.write('	\hyperref[sec:nm_sec_%s]{%s-}\hyperref[sec:nm_sec_%s]{%s}& %s \\\\\n'%(opt_name, opt_name[0:idx].replace('_','\_'), opt_name, opt_name[idx:].replace('_','\_'), opt_description))
			else:
				latex.write('	\hyperref[sec:nm_sec_%s]{%s} & %s \\\\\n'%(opt_name, opt_name.replace('_','\_'), opt_description))
			latex.write('	\hline\n')

		latex.write('\end{longtable}\n')
		latex.write('\end{center}\n')
		latex.write('}\n')
#}}}

def write_variable_table(latex, registry, mode):#{{{
	latex_missing_string = '{\\bf \color{red} MISSING}'
	variable_table_header = '{\\bf Name} & {\\bf Description}'

	latex.write('\section[Variable definitions]{\hyperref[chap:variable_sections]{Variable definitions}}\n')
	latex.write('\label{sec:%s_variable_tables}\n'%mode)
	latex.write('Embedded links point to more detailed variable information in the appendix.\n')
	for var_struct in registry.iter("var_struct"):
		struct_name = var_struct.attrib['name']
		latex.write('\subsection[%s]{\hyperref[sec:var_sec_%s]{%s}}\n'%(struct_name.replace('_','\_'),struct_name,struct_name.replace('_','\_')))
		latex.write('\label{subsec:%s_var_tab_%s}\n'%(mode, struct_name))

		try:
			junk_file = open('%s/%s_struct.tex'%(options.latex_dir,struct_name), 'r')
			latex.write('\input{%s/%s_struct.tex}\n'%(options.latex_path, struct_name))
			junk_file.close()
		except:
			latex.write('')

		latex.write('\\vspace{0.5in}\n')
		latex.write('{\small\n')
		latex.write('\\begin{center}\n')
		latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
		latex.write('	\hline\n')
		latex.write('	%s \\endfirsthead\n'%variable_table_header)
		latex.write('	\hline \n')
		latex.write('	%s (Continued) \\endhead\n'%variable_table_header)
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

					idx = break_string(var_name)
					if idx > -1:
						latex.write('	\hyperref[subsec:var_sec_%s_%s]{%s-}\hyperref[subsec:var_sec_%s_%s]{%s}  & %s \\\\\n'%(struct_name, var_name, var_name[0:idx].replace('_','\_'), struct_name, var_name, var_name[idx:].replace('_','\_'), var_description))
					else:
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

				idx = break_string(var_name)
				if idx > -1:
					latex.write('	\hyperref[subsec:var_sec_%s_%s]{%s-}\hyperref[subsec:var_sec_%s_%s]{%s  }& %s \\\\\n'%(struct_name, var_name, var_name[0:idx].replace('_','\_'), struct_name, var_name, var_name[idx:].replace('_','\_'), var_description))
				else:
					latex.write('	\hyperref[subsec:var_sec_%s_%s]{%s} & %s \\\\\n'%(struct_name, var_name, var_name.replace('_','\_'), var_description))
				latex.write('	\hline\n')

		latex.write('\end{longtable}\n')
		latex.write('\end{center}\n')
		latex.write('}\n')
#}}}

def write_namelist_sections(latex, sorted_opts, forward_registry, analysis_registry):#{{{
	latex_missing_string = '{\\bf \color{red} MISSING}'

	#latex.write('\chapter[Namelist options]{\hyperref[chap:namelist_tables]{Namelist options}}\n')
	latex.write('\chapter[Namelist options]{Namelist options}\n')
	latex.write('\label{chap:namelist_sections}\n')
#	latex.write('Embedded links point to information in chapter \\ref{chap:namelist_tables}\n')

	for opt in sorted_opts:
		found = False
		in_forward = False
		in_analysis = False
		forward_rec_name = ""
		analysis_rec_name = ""

		# Search forward registry
		for nml_rec in forward_registry.iter("nml_record"):#{{{
			for nml_opt in nml_rec.iter("nml_option"):
				opt_name = nml_opt.attrib["name"]

				if(in_forward == False and opt_name == opt):
					in_forward = True
					forward_rec_name = nml_rec.attrib["name"]
					if(not found):
						found = True
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
#}}}

		# Search analysis registry if not found yet
		for nml_rec in analysis_registry.iter("nml_record"):#{{{
			for nml_opt in nml_rec.iter("nml_option"):
				opt_name = nml_opt.attrib["name"]

				if(in_analysis == False and opt_name == opt):
					in_analysis = True
					analysis_rec_name = nml_rec.attrib["name"]
					if(not found):
						found = True
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
#}}}

		# If option has been found, write it out as a section.
		if(found):
			opt_name = opt
			#latex.write('\section[%s]{\hyperref[sec:nm_tab_%s]{%s}}\n'%(opt_name.replace('_','\_'),rec_name,opt_name.replace('_','\_')))
			latex.write('\section[%s]{%s}\n'%(opt_name.replace('_','\_'),opt_name.replace('_','\_')))
			latex.write('\label{sec:nm_sec_%s}\n'%opt_name)
			latex.write('\\begin{center}\n')
			latex.write('\\begin{longtable}{| p{2.0in} || p{4.0in} |}\n')
			latex.write('    \hline\n')
			latex.write('    In build modes: & ')
			if(in_forward):
				latex.write('\hyperref[subsec:forward_nm_tab_%s]{forward} '%(forward_rec_name))
			if(in_analysis):
				latex.write('\hyperref[subsec:analysis_nm_tab_%s]{analysis} '%(analysis_rec_name))
			latex.write('\\\\\n')

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
#}}}

def write_variable_sections(latex, sorted_structs, forward_registry, analysis_registry):#{{{
	latex_missing_string = '{\\bf \color{red} MISSING}'

	#latex.write('\chapter[Variable definitions]{\hyperref[chap:variable_tables]{Variable definitions}}\n')
	latex.write('\chapter[Variable definitions]{Variable definitions}\n')
	latex.write('\label{chap:variable_sections}\n')
#	latex.write('Embedded links point to information in chapter \\ref{chap:variable_tables}\n')

	for struct in sorted_structs:
		##latex.write('\section[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(struct.replace('_','\_'),struct, struct.replace('_','\_')))
		latex.write('\section[%s]{%s}\n'%(struct.replace('_','\_'), struct.replace('_','\_')))
		latex.write('\label{sec:var_sec_%s}\n'%struct)

		unique_vars = [];
		# Determine all variables in the current var struct from the forward mode
		for var_struct in forward_registry.iter("var_struct"):#{{{
			struct_name = var_struct.attrib["name"]
			if(struct_name == struct):
				for var_arr in var_struct.iter("var_array"):
					for var in var_arr.iter("var"):
						name = var.attrib["name"]
						if(unique_vars.count(name) == 0):
							unique_vars.append(name)
				for var in var_struct.iter("var"):
					name = var.attrib["name"]
					if(unique_vars.count(name) == 0):
						unique_vars.append(name)
#}}}

		# Determine all variables in the current var struct from the analysis mode
		for var_struct in analysis_registry.iter("var_struct"):#{{{
			struct_name = var_struct.attrib["name"]
			if(struct_name == struct):
				for var_arr in var_struct.iter("var_array"):
					for var in var_arr.iter("var"):
						name = var.attrib["name"]
						if(unique_vars.count(name) == 0):
							unique_vars.append(name)
				for var in var_struct.iter("var"):
					name = var.attrib["name"]
					if(unique_vars.count(name) == 0):
						unique_vars.append(name)
			#}}}

		sorted_vars = sorted(unique_vars)
		del unique_vars

		for var_name in sorted_vars:
			found = False
			in_forward = False
			in_analysis = False
			in_var_array = False
			
			# Try to extract var from forward mode
			for var_struct in forward_registry.iter("var_struct"): #{{{
				struct_name = var_struct.attrib["name"]
				if(struct_name == struct):
					for var_arr in var_struct.iter("var_array"):
						for var in var_arr.iter("var"):
							name = var.attrib["name"]
							if(name == var_name):
								in_forward = True

								if(not found):
									found = True
									in_var_array = True

									struct_time_levs = var_struct.attrib['time_levs']
									var_arr_name = var_arr.attrib["name"]
									var_type = var_arr.attrib["type"]
									var_dims = var_arr.attrib["dimensions"]
									try:
									    var_time_levels = var_arr.attrib["time_levs"]
									except:
									    var_time_levels = var_struct.attrib["time_levs"]

									# Extract var persistence#{{{
									try:
										var_persistence = var_arr.attrib['persistence']
									except:
										var_persistence = 'persistent'
#}}}

									# Extract var units#{{{
									try:
										var_units = var.attrib['units']
										if var_units == "":
											var_units = latex_missing_string
										else:
											var_units = "$%s$"%var_units.replace(' ', '$ $')
									except:
										var_units = latex_missing_string
#}}}

									var_arr_group = var.attrib["array_group"]

									# Extract name in code, and build var_path#{{{
									try:
										var_name_in_code = var.attrib['name_in_code']
									except:
										var_name_in_code = var_name

									var_path = "%s"%(var_arr_name)
#}}}

									# Extract var description#{{{
									try:
										var_description = var.attrib['description']
									except:
										var_description = latex_missing_string.replace('_','\_')

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
#}}}

					for var in var_struct.iter("var"):
						name = var.attrib["name"]
						if(name == var_name):
							in_forward = True

							if(not found):
								found = True
								in_var_array = False
								struct_time_levs = var_struct.attrib['time_levs']
								var_type = var.attrib["type"]
								var_dims = var.attrib["dimensions"]
								try:
								    var_time_levels = var.attrib["time_levs"]
								except:
								    var_time_levels = var_struct.attrib["time_levs"]

								# Extract var persistence#{{{
								try:
									var_persistence = var_arr.attrib['persistence']
								except:
									var_persistence = 'persistent'
#}}}

								# Extract var units#{{{
								try:
									var_units = var.attrib['units']
									if var_units == "":
										var_units = latex_missing_string
									else:
										var_units = "$%s$"%var_units.replace(' ', '$ $')
								except:
									var_units = latex_missing_string
#}}}

								# Extract name in code, and build var_path#{{{
								try:
									var_name_in_code = var.attrib['name_in_code']
								except:
									var_name_in_code = var_name

								var_path = "%s"%(var_name)
#}}}

								# Extract var description#{{{
								try:
									var_description = var.attrib['description']
								except:
									var_description = latex_missing_string.replace('_','\_')

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
#}}}
#}}}

			# Try to extract var from analysis mode
			for var_struct in analysis_registry.iter("var_struct"): #{{{
				struct_name = var_struct.attrib["name"]
				if(struct_name == struct):
					for var_arr in var_struct.iter("var_array"):
						for var in var_arr.iter("var"):
							name = var.attrib["name"]
							if(name == var_name):
								in_analysis = True

								if(not found):
									found = True
									in_var_array = True

									struct_time_levs = var_struct.attrib['time_levs']
									var_arr_name = var_arr.attrib["name"]
									var_type = var_arr.attrib["type"]
									var_dims = var_arr.attrib["dimensions"]
									try:
									    var_time_levels = var_arr.attrib["time_levs"]
									except:
									    var_time_levels = var_struct.attrib["time_levs"]

									# Extract var persistence#{{{
									try:
										var_persistence = var_arr.attrib['persistence']
									except:
										var_persistence = 'persistent'
#}}}

									# Extract var units#{{{
									try:
										var_units = var.attrib['units']
										if var_units == "":
											var_units = latex_missing_string
										else:
											var_units = "$%s$"%var_units.replace(' ', '$ $')
									except:
										var_units = latex_missing_string
#}}}

									var_arr_group = var.attrib["array_group"]

									# Extract name in code, and build var_path#{{{
									try:
										var_name_in_code = var.attrib['name_in_code']
									except:
										var_name_in_code = var_name

									if int(struct_time_levs) > 1:
										var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code.replace('_','\_'))
										var_path = "%s"%(var_arr_name)
									else:
										var_index = "domain %% blocklist %% %s %% index_%s"%(struct_name, var_name_in_code.replace('_','\_'))
										var_path = "%s"%(var_arr_name)
#}}}

									# Extract var description#{{{
									try:
										var_description = var.attrib['description']
									except:
										var_description = latex_missing_string.replace('_','\_')

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
#}}}

					for var in var_struct.iter("var"):
						name = var.attrib["name"]
						if(name == var_name):
							in_analysis = True

							if(not found):
								found = True
								in_var_array = False
								struct_time_levs = var_struct.attrib['time_levs']
								var_type = var.attrib["type"]
								var_dims = var.attrib["dimensions"]
								try:
								    var_time_levels = var.attrib["time_levs"]
								except:
								    var_time_levels = var_struct.attrib["time_levs"]

								# Extract var persistence#{{{
								try:
									var_persistence = var_arr.attrib['persistence']
								except:
									var_persistence = 'persistent'
#}}}

								# Extract var units#{{{
								try:
									var_units = var.attrib['units']
									if var_units == "":
										var_units = latex_missing_string
									else:
										var_units = "$%s$"%var_units.replace(' ', '$ $')
								except:
									var_units = latex_missing_string
#}}}

								# Extract name in code, and build var_path#{{{
								try:
									var_name_in_code = var.attrib['name_in_code']
								except:
									var_name_in_code = var_name

								if int(struct_time_levs) > 1:
									var_path = "domain %% blocklist %% %s %% time_levs(:) %% %s %% %s"%(struct_name, struct_name, var_name)
								else:
									var_path = "domain %% blocklist %% %s %% %s"%(struct_name, var_name)
#}}}

								# Extract var description#{{{
								try:
									var_description = var.attrib['description']
								except:
									var_description = latex_missing_string.replace('_','\_')

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
#}}}
#}}}

			# Build stream list from forward mode
			forward_streams = ""#{{{
			if(in_forward):
				for streams in forward_registry.iter("streams"):
					for stream in streams.iter("stream"):
						for var in stream.iter("var"):
							name = var.attrib["name"]
							if(name == var_name):
								forward_streams = "%s %s"%(forward_streams, stream.attrib["name"])
								#}}}

			# Build stream list from analysis mode
			analysis_streams = ""#{{{
			if(in_analysis):
				for streams in analysis_registry.iter("streams"):
					for stream in streams.iter("stream"):
						for var in stream.iter("var"):
							name = var.attrib["name"]
							if(name == var_name):
								analysis_streams = "%s %s"%(analysis_streams, stream.attrib["name"])
								#}}}

			if(found):
				struct_name = struct
				#latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n'%(var_name.replace('_','\_'),struct_name, var_name.replace('_','\_')))
				latex.write('\subsection[%s]{%s}\n'%(var_name.replace('_','\_'), var_name.replace('_','\_')))
				latex.write('\label{subsec:var_sec_%s_%s}\n'%(struct_name,var_name))
				# Tabular Format:
				latex.write('\\begin{center}\n')
				latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
				latex.write('        \hline \n')
				latex.write('        In build modes: & ')
				if(in_forward):
					latex.write('\hyperref[subsec:forward_var_tab_%s]{forward} '%(struct))
				if(in_analysis):
					latex.write('\hyperref[subsec:analysis_var_tab_%s]{analysis} '%(struct))
				latex.write('\\\\\n')
				latex.write('        \hline \n')
				latex.write('        Type: & %s \\\\\n'%var_type)
				latex.write('        \hline \n')
				latex.write('        Units: & %s \\\\\n'%var_units)
				latex.write('        \hline \n')
				latex.write('        Dimension: & %s \\\\\n'%var_dims)
				latex.write('        \hline \n')
				latex.write('        Persistence: & %s \\\\\n'%var_persistence)
				latex.write('        \hline \n')
				latex.write('        Number of time levels: & %s \\\\\n'%var_time_levels)
				latex.write('        \hline \n')

				if(in_var_array):
					latex.write("		 Index in `%s' Array: & `index\_%s' in `%s' pool \\\\\n"%(var_path.replace('_','\_'), var_name.replace('_', '\_'), struct_name.replace('_','\_')))
					latex.write('		 \hline \n')
				pool_path="`%s' in `%s' pool"%(var_path, struct_name)
				latex.write('            Pool path: & %s \\\\\n'%pool_path.replace('_', '\_'))
				latex.write('		 \hline \n')
				if(in_var_array):
					latex.write('		 Array Group: & %s \\\\\n'%var_arr_group.replace('_','\_'))
					latex.write('		 \hline \n')
				latex.write('    \caption{%s: %s}\n'%(var_name.replace('_','\_'),var_description))
				latex.write('\end{longtable}\n')
				latex.write('\end{center}\n')

#}}}

if not options.forward_registry_path:
	parser.error("The forward mode's Registry file is required")

if not options.analysis_registry_path:
	parser.error("The analysis mode's Registry file is required")

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

forward_registry_path = options.forward_registry_path
analysis_registry_path = options.analysis_registry_path

try:
	forward_registry_tree = ET.parse(forward_registry_path)
except:
	parser.error('%s does not exist or is not parsable. Exiting.'%forward_registry_path)

try:
	analysis_registry_tree = ET.parse(analysis_registry_path)
except:
	parser.error('%s does not exist or is not parsable. Exiting.'%analysis_registry_path)

forward_registry = forward_registry_tree.getroot()
analysis_registry = analysis_registry_tree.getroot()

# Build dimension lists
forward_dims  = [];
analysis_dims = [];
unique_dims   = [];

for dims in forward_registry.iter("dims"):
	for dim in dims.iter("dim"):
		dim_name = dim.attrib["name"]
		forward_dims.append(dim_name)
		if(unique_dims.count(dim_name) == 0):
			unique_dims.append(dim_name)

for dims in analysis_registry.iter("dims"):
	for dim in dims.iter("dim"):
		dim_name = dim.attrib["name"]
		analysis_dims.append(dim_name)
		if(unique_dims.count(dim_name) == 0):
			unique_dims.append(dim_name)

sorted_dims = sorted(unique_dims)
del unique_dims

# Build structure lists
forward_structs = [];
analysis_structs = [];
unique_structs = [];

for struct in forward_registry.iter("var_struct"):
	struct_name = struct.attrib["name"]
	forward_structs.append(struct_name)
	if(unique_structs.count(struct_name) == 0):
		unique_structs.append(struct_name)

for struct in analysis_registry.iter("var_struct"):
	struct_name = struct.attrib["name"]
	analysis_structs.append(struct_name)
	if(unique_structs.count(struct_name) == 0):
		unique_structs.append(struct_name)

sorted_structs = sorted(unique_structs)
del unique_structs


# Build Namelist lists
forward_opts = [];
analysis_opts = [];
unique_opts = [];

for nml_rec in forward_registry.iter("nml_record"):
	for nml_opt in nml_rec.iter("nml_option"):
		opt_name = nml_opt.attrib["name"]
		forward_opts.append(opt_name)
		if(unique_opts.count(opt_name) == 0):
			unique_opts.append(opt_name)

for nml_rec in analysis_registry.iter("nml_record"):
	for nml_opt in nml_rec.iter("nml_option"):
		opt_name = nml_opt.attrib["name"]
		analysis_opts.append(opt_name)
		if(unique_opts.count(opt_name) == 0):
			unique_opts.append(opt_name)

sorted_opts = sorted(unique_opts)
del unique_opts

# Write file that defines version string for model.
latex = open('define_version.tex', 'w+')
try:
	version_string = forward_registry.attrib['version']
except:
	version_string = '{\\bf MISSING}'
latex.write('\\newcommand{\\version}{%s}\n'%version_string)
latex.close()

# Write file to include for forward mode
# It should have sections for dimensions, namelist options/records, and
# variable and their structures.
latex = open('mode_forward_sections.tex', 'w+')
mode = 'forward'
write_dimension_table(latex, forward_registry, mode)
write_namelist_table(latex, forward_registry, mode)
write_variable_table(latex, forward_registry, mode)
latex.close()

# Write file to include for analysis mode
# It should have sections for dimensions, namelist options/records, and
# variable and their structures.
latex = open('mode_analysis_sections.tex', 'w+')
mode = 'analysis'
write_dimension_table(latex, analysis_registry, mode)
write_namelist_table(latex, analysis_registry, mode)
write_variable_table(latex, analysis_registry, mode)
latex.close()


latex = open('namelist_sections.tex', 'w+')
write_namelist_sections(latex, sorted_opts, forward_registry, analysis_registry)
latex.close()

latex = open('variable_sections.tex', 'w+')
write_variable_sections(latex, sorted_structs, forward_registry, analysis_registry)
latex.close()


