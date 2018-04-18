#!/usr/bin/env python

"""
This script parses a MPAS Registry.xml file to generates documentation for a
users or developers guide.

Typical usage is as follows::

    # set the core, one of ocean, landice, cice, etc.
    export CORE=<core>
    # Set your repo directories:
    export MPAS_REPO=~/repos/MPAS
    export MPAS_TOOLS_REPO=~/repos/MPAS-Tools
    export MPAS_DOCUMENTS_REPO=~/repos/MPAS-Documents
    cd $MPAS_REPO
    # Compile MPAS so you have a src/core_ocean/Registry_processed.xml file.
    # Change the compiler as needed.
    make CORE=$CORE gfortran
    cd $MPAS_DOCUMENTS_REPO/users_guide/$CORE
    # clean up blank lines at the top of the XML file
    sed '/./,$!d' $MPAS_REPO/src/core_${CORE}/Registry_processed.xml > \
      Registry_cleaned.xml
    $MPAS_TOOLS_REPO/python_scripts/namelist_generation/parse_xml_registry.py \
      -f Registry_cleaned.xml -d section_descriptions \
      -p ${CORE}/section_descriptions
    cd ..
    make clean CORE=$CORE
    make CORE=$CORE

The -f flag points to the processed registry file (typically with a full path).

The -d flag points to the local or full path to .tex files that containing
section descriptions for providing additional information in the output latex
documentation.

Section descriptions are required to be named whatever the section is. For
example, in a namelist, there might be a namelist record named
"&time_management". The script searches the directory listed with the -d
flag for a latex file named time_management.tex, and adds an input line to
the output latex documentation to include this file.

The -p flag specifies the relative path inside the latex documentation where
the file should be input from. As an example, one might
run it as follows to generate the ocean core's documentation::

    ./parse_xml_registry.xml -f mpas_root/src/core_ocean/Registry.xml \
        -d mpas_doc_root/users_guide/ocean/section_descriptions \
        -p ocean/section_descriptions

On output, several files are created which are listed below.
    namelist.input.generated - A default namelist.input file for the core that
                               owns the Registry.xml file.
    dimensions.tex - A tabulated description of the dimensions for the core.
    namelist_table_documentation.tex - A tabulated description of the namelist
                                       options for the core.
    namelist_section_documentation.tex - A more detailed section format
                                         description of the namelist options
                                         for the core.
    variable_table_documentation.tex - A tabulated description of the variables
                                       in the core.
    variable_section_documentation.tex - A more detailed section formate
                                         description of the variable in the
                                         core.
    define_version.tex - A simple file which can be included to define \version
                         inside the users guide.

Authors:
========
Doug Jacobsen, Xylar Asay-Davis
"""


import os
from optparse import OptionParser
import xml.etree.ElementTree as ET
from collections import OrderedDict
from PIL import ImageFont
import pkg_resources


def break_string(string, maxLength=150., font='cmunrm.otf', fontSize=10):
    # {{{

    # Note: max_length is in points, so 144. corresponds to 2 inches, the
    # column width for namelist and variable names in tables in the user guide

    # font defaults to LaTex font (Computer Modern), and user guide font size
    # in tables

    # if an absolute path to the font was not supplied, look relative to this
    # script
    if not os.path.isabs(font):
        font = pkg_resources.resource_filename(__name__, font)

    font = ImageFont.truetype(font, fontSize)
    size = font.getsize(string)
    if size[0] <= maxLength:
        # no need to split
        return None

    bestBreakPoints = []

    # first alpha-numeric character after a non-alpha-numeric character
    for index in range(1, len(string)):
        if not string[index-1].isalnum() and string[index].isalnum():
            bestBreakPoints.append(index)

    # find uppercase following lowercase or number
    for index in range(1, len(string)):
        if string[index-1].isalnum() and string[index-1].islower() \
                and string[index].isalpha() and string[index].isupper():
            bestBreakPoints.append(index)

    bestBreakPoints.append(len(string))

    bestBreakPoints = sorted(bestBreakPoints)

    for index in range(1, len(bestBreakPoints)):
        breakPoint = bestBreakPoints[index]
        size = font.getsize(string[:breakPoint])
        if size[0] > maxLength:
            breakPoint = bestBreakPoints[index-1]
            return breakPoint

    # there is no good break point so we have to find an arbitrary one
    print "Warning: no good breakpoint found for {}".format(string)
    for breakPoint in range(1, len(string)+1):
        breakPoint = bestBreakPoints[index]
        size = font.getsize(string[:breakPoint])
        if size[0] > maxLength:
            breakPoint = breakPoint-1
            return breakPoint

    raise ValueError("Could not find a breakpoint for {}".format(string))
    # }}}


def write_namelist_input_generated():
    # Write default namelist
    namelist = open('namelist.input.generated', 'w')
    for nml_rec in registry.iter("nml_record"):
        namelist.write('&%s\n' % nml_rec.attrib['name'])
        for nml_opt in nml_rec.iter("nml_option"):
            if nml_opt.attrib['type'] == "character":
                namelist.write('\t%s = "%s"\n' % (
                        nml_opt.attrib['name'],
                        nml_opt.attrib['default_value']))
            else:
                namelist.write('\t%s = %s\n' % (
                        nml_opt.attrib['name'],
                        nml_opt.attrib['default_value']))

        namelist.write('/\n')


def escape_underscore(string):
    has_math_mode = (string.find('$') != -1)
    if has_math_mode:
        dim_desc_split = string.split("$")
        replace = True
        string = ""
        for part in dim_desc_split:
            if replace:
                part = part.replace('_', '\_')
                string = "{}{}".format(string, part)
                replace = False
            else:
                string = "{}${}$".format(string, part)
                replace = True
    else:
        string = string.replace('_', '\_')
    return string


def get_attrib(element, attributeName, missingValue=None):
    if missingValue is None:
        missingValue = latex_missing_string
    try:
        attrib = element.attrib[attributeName]
    except KeyError:
        attrib = missingValue
    if attrib == "":
        attrib = missingValue
    return attrib


def get_units(element):
    units = get_attrib(element, 'units')
    if units != latex_missing_string:
        # units with the siunitx package
        units = "\si{{{}}}".format(units.replace(' ', '.'))
    units = escape_underscore(units)
    return units


def get_description(element):
    description = get_attrib(element, 'description')
    description = escape_underscore(description)
    return description


def get_linked_name(name, link):
    indices = []
    index = 0
    while True:
        newIndex = break_string(name[index:])
        if newIndex is None:
            break
        index += newIndex
        indices.append(index)

    indices.append(len(name))
    newName = escape_underscore(name[0:indices[0]])
    for start, end in zip(indices[0:-1], indices[1:]):
        namePiece = escape_underscore(name[start:end])
        newName = '{}\\-{}'.format(newName, namePiece)

    return '\hyperref[subsec:%s]{%s}' % (link, newName)


def write_var_struct_to_table(latex, var_struct, struct_name):
    for node in var_struct:
        if node.tag == 'var_struct':
            write_var_struct_to_table(latex, node, struct_name)
        elif node.tag == 'var_array':
            write_var_array_to_table(latex, node, struct_name)
        elif node.tag == 'var':
            write_var_to_table(latex, node, struct_name)


def write_var_array_to_table(latex, var_array, struct_name):
    for var in var_array.iter("var"):
        write_var_to_table(latex, var, struct_name)


def write_var_to_table(latex, var, struct_name):
    var_name = var.attrib['name']
    var_description = get_description(var)

    link = 'var_sec_{}_{}'.format(struct_name, var_name)
    linkedName = get_linked_name(var_name, link)

    latex.write('    {} & {} \\\\\n'.format(linkedName,
                                            var_description))
    latex.write('    \hline\n')


def get_var_structs():
    # use a dictionary to create lists of all top-level var_structs with the
    # same name (e.g. state, tracers, mesh)
    var_structs = OrderedDict()
    for var_struct in registry:
        if var_struct.tag != "var_struct":
            continue
        struct_name = var_struct.attrib['name']
        if struct_name in var_structs.keys():
            var_structs[struct_name].append(var_struct)
        else:
            var_structs[struct_name] = [var_struct]
    return var_structs


def write_var_struct_section(latex, var_struct, struct_name, has_time):
    for node in var_struct:
        if node.tag == 'var_struct':
            write_var_struct_section(latex, node, struct_name, has_time)
        elif node.tag == 'var_array':
            write_var_array_section(latex, node, struct_name, has_time)
        elif node.tag == 'var':
            write_var_section(latex, node, struct_name, has_time)


def write_var_array_section(latex, var_array, struct_name, has_time):
    for var in var_array.iter("var"):
        write_var_section(latex, var, struct_name, has_time, var_array)


def write_var_section(latex, var, struct_name, has_time, var_array=None):
    var_name = var.attrib['name']
    var_name_escaped = escape_underscore(var_name)
    if var_array is None:
        var_type = var.attrib['type']
        dimensions = var.attrib['dimensions']
    else:
        var_arr_name = escape_underscore(var_array.attrib['name'])
        var_type = var_array.attrib['type']
        dimensions = var_array.attrib['dimensions']

    persistence = get_attrib(var, "persistence", missingValue='persistent')
    name_in_code = get_attrib(var, "name_in_code", missingValue=var_name)
    units = get_units(var)
    description = get_description(var)

    if has_time:
        var_path = "domain % blocklist % {} % time_levs(:) % {} % {}".format(
                struct_name, struct_name, var_name)
    else:
        var_path = "domain % blocklist % {} % {}".format(struct_name, var_name)

    var_path = escape_underscore(var_path).replace('%', '\%')

    latex.write('\subsection[%s]{\hyperref[sec:var_tab_%s]{%s}}\n' % (
            var_name_escaped, struct_name, var_name_escaped))
    latex.write('\label{subsec:var_sec_%s_%s}\n' % (struct_name, var_name))
    # Tabular Format:
    latex.write('\\begin{center}\n')
    latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
    latex.write('        \hline \n')
    latex.write('        Type: & %s \\\\\n' % var_type)
    latex.write('        \hline \n')
    latex.write('        Units: & %s \\\\\n' % units)
    latex.write('        \hline \n')
    latex.write('        Dimension: & %s \\\\\n' % dimensions)
    latex.write('        \hline \n')
    latex.write('        Persistence: & %s \\\\\n' % (persistence))
    latex.write('        \hline \n')

    if var_array is not None:
        array_group = escape_underscore(var.attrib['array_group'])
        index = "domain % blocklist % {} % index_{}".format(struct_name,
                                                            name_in_code)
        index = escape_underscore(index).replace('%', '\%')

        latex.write('         Index in %s Array: & %s \\\\\n' % (var_arr_name,
                                                                 index))
        latex.write('         \hline \n')

    latex.write('         Location in code: & %s \\\\\n' % (var_path))
    latex.write('         \hline \n')

    if var_array is not None:
        latex.write('         Array Group: & %s \\\\\n' % (array_group))
        latex.write('         \hline \n')

    latex.write('    \caption{%s: %s}\n' % (var_name_escaped, description))
    latex.write('\end{longtable}\n')
    latex.write('\end{center}\n')


def write_dimension_table_documentation():
    # Write dimension table documentation latex file.
    latex = open('dimension_table_documentation.tex', 'w')
    latex.write('\chapter{Dimensions}\n')
    latex.write('\label{chap:dimensions}\n')
    latex.write('{\small\n')
    latex.write('\\begin{center}\n')
    latex.write('\\begin{longtable}{| p{1.0in} || p{1.0in} | p{4.0in} |}\n')
    latex.write('    \hline \n')
    latex.write('    {} \\endfirsthead\n'.format(dimension_table_header))
    latex.write('    \hline \n')
    latex.write('    {} (Continued) \\endhead\n'.format(
            dimension_table_header))
    latex.write('    \hline \n')
    latex.write('    \hline \n')
    for dims in registry.iter("dims"):
        for dim in dims.iter("dim"):
            name = dim.attrib['name']
            name = escape_underscore(name)
            units = get_units(dim)
            description = get_description(dim)

            latex.write('    {} & {} & {} \\\\ \n'.format(
                    name, units, description))
            latex.write('    \hline\n')

    latex.write('\end{longtable}\n')
    latex.write('\end{center}\n')
    latex.write('}\n')
    latex.close()


def write_namelist_table_documentation():
    # Write namelist table documentation latex file.
    latex = open('namelist_table_documentation.tex', 'w')
    latex.write('\chapter[Namelist options]{\hyperref[chap:namelist_sections]'
                '{Namelist options}}\n')
    latex.write('\label{chap:namelist_tables}\n')
    latex.write('Embedded links point to more detailed namelist information '
                'in the appendix.\n')
    for nml_rec in registry.iter("nml_record"):
        rec_name = nml_rec.attrib['name']
        rec_name_escaped = escape_underscore(rec_name)
        latex.write('\section[%s]{\hyperref[sec:nm_sec_%s]{%s}}\n' % (
                rec_name_escaped, rec_name, rec_name_escaped))
        latex.write('\label{sec:nm_tab_%s}\n' % (rec_name))

        # Add input line if file exists.
        if os.path.exists('%s/%s.tex' % (options.latex_dir, rec_name)):
            latex.write('\input{%s/%s.tex}\n' % (options.latex_path, rec_name))
        else:
            print 'Warning, namelist description latex file not found: ' \
                  '%s/%s.tex' % (options.latex_dir, rec_name)
            latex.write('')

        latex.write('\\vspace{0.5in}\n')
        latex.write('{\small\n')
        latex.write('\\begin{center}\n')
        latex.write('\\begin{longtable}{| p{2.0in} || p{4.0in} |}\n')
        latex.write('    \hline\n')
        latex.write('    %s \\endfirsthead\n' % namelist_table_header)
        latex.write('    \hline \n')
        latex.write('    %s (Continued) \\endhead\n' % namelist_table_header)
        latex.write('    \hline\n')
        latex.write('    \hline\n')

        for nml_opt in nml_rec.iter("nml_option"):
            name = nml_opt.attrib['name']

            description = get_description(nml_opt)

            link = 'nm_sec_{}'.format(name)
            linkedName = get_linked_name(name, link)

            latex.write('    {} & {} \\\\\n'.format(linkedName, description))
            latex.write('    \hline\n')

        latex.write('\end{longtable}\n')
        latex.write('\end{center}\n')
        latex.write('}\n')
    latex.close()


def write_namelist_section_documentation():
    # Write namelist section documentation latex file.
    latex = open('namelist_section_documentation.tex', 'w')
    latex.write('\chapter[Namelist options]{\hyperref[chap:namelist_tables]'
                '{Namelist options}}\n')
    latex.write('\label{chap:namelist_sections}\n')
    latex.write('Embedded links point to information in chapter '
                '\\ref{chap:namelist_tables}\n')
    for nml_rec in registry.iter("nml_record"):
        rec_name = nml_rec.attrib["name"]
        rec_name_escaped = escape_underscore(rec_name)

        latex.write('\section[%s]{\hyperref[sec:nm_tab_%s]{%s}}\n' % (
                rec_name_escaped, rec_name, rec_name_escaped))
        latex.write('\label{sec:nm_sec_%s}\n' % rec_name)

        for nml_opt in nml_rec.iter("nml_option"):
            name = nml_opt.attrib["name"]
            name_escaped = escape_underscore(name)
            opt_type = escape_underscore(nml_opt.attrib["type"])
            default_value = escape_underscore(get_attrib(nml_opt,
                                                         "default_value"))
            possible_values = escape_underscore(get_attrib(nml_opt,
                                                           "possible_values"))
            units = get_units(nml_opt)
            description = get_description(nml_opt)

            try:
                opt_icepack_name = nml_opt.attrib["icepack_name"]
            except KeyError:
                opt_icepack_name = None

            latex.write('\subsection[%s]{\hyperref[sec:nm_tab_%s]{%s}}\n' % (
                    name_escaped, rec_name, name_escaped))
            latex.write('\label{subsec:nm_sec_%s}\n' % name)
            latex.write('\\begin{center}\n')
            latex.write('\\begin{longtable}{| p{2.0in} || p{4.0in} |}\n')
            latex.write('    \hline\n')
            latex.write('    Type: & %s \\\\\n' % opt_type)
            latex.write('    \hline\n')
            latex.write('    Units: & %s \\\\\n' % units)
            latex.write('    \hline\n')
            latex.write('    Default Value: & %s \\\\\n' % default_value)
            latex.write('    \hline\n')
            latex.write('    Possible Values: & %s \\\\\n' % possible_values)
            latex.write('    \hline\n')
            if (opt_icepack_name is not None):
                latex.write('    Icepack name: & \\verb+%s+ \\\\\n' %
                            opt_icepack_name)
                latex.write('    \hline\n')
            latex.write('    \caption{%s: %s}\n' % (name_escaped, description))
            latex.write('\end{longtable}\n')
            latex.write('\end{center}\n')
    latex.close()


def write_variable_table_documentation():

    # Write variable table documentation latex file
    latex = open('variable_table_documentation.tex', 'w')
    latex.write('\chapter[Variable definitions]'
                '{\hyperref[chap:variable_sections]'
                '{Variable definitions}}\n')
    latex.write('\label{chap:variable_tables}\n')
    latex.write('Embedded links point to more detailed variable information '
                'in the appendix.\n')

    var_structs = get_var_structs()

    for struct_name, var_struct_list in var_structs.items():
        struct_name_escaped = escape_underscore(struct_name)
        latex.write('\section[%s]{\hyperref[sec:var_sec_%s]{%s}}\n' % (
                struct_name_escaped, struct_name, struct_name_escaped))
        latex.write('\label{sec:var_tab_%s}\n' % struct_name)

        if os.path.exists('%s/%s_struct.tex' % (options.latex_dir,
                                                struct_name)):
            latex.write('\input{%s/%s_struct.tex}\n' % (options.latex_path,
                                                        struct_name))
        else:
            print 'Warning, variable section description latex file not ' \
                'found:  %s/%s_struct.tex' % (options.latex_dir, struct_name)
            latex.write('')

        latex.write('\\vspace{0.5in}\n')
        latex.write('{\small\n')
        latex.write('\\begin{center}\n')
        latex.write('\\begin{longtable}{| p{2.0in} | p{4.0in} |}\n')
        latex.write('    \hline\n')
        latex.write('    %s \\endfirsthead\n' % variable_table_header)
        latex.write('    \hline \n')
        latex.write('    %s (Continued) \\endhead\n' % variable_table_header)
        latex.write('    \hline\n')

        for var_struct in var_struct_list:
            write_var_struct_to_table(latex, var_struct, struct_name)

        latex.write('\end{longtable}\n')
        latex.write('\end{center}\n')
        latex.write('}\n')
    latex.close()


def write_variable_section_documentation():

    # Write variable section documentation latex file
    latex = open('variable_section_documentation.tex', 'w')
    latex.write('\chapter[Variable definitions]'
                '{\hyperref[chap:variable_tables]'
                '{Variable definitions}}\n')
    latex.write('\label{chap:variable_sections}\n')
    latex.write('Embedded links point to information in chapter '
                '\\ref{chap:variable_tables}\n')

    var_structs = get_var_structs()

    for struct_name, var_struct_list in var_structs.items():
        struct_name_escaped = escape_underscore(struct_name)

        latex.write('\section[%s]{\hyperref[sec:var_tab_%s]{%s}}\n' % (
                struct_name_escaped, struct_name, struct_name_escaped))
        latex.write('\label{sec:var_sec_%s}\n' % struct_name)

        for var_struct in var_struct_list:
            try:
                struct_time_levs = var_struct.attrib['time_levs']
                has_time = int(struct_time_levs) > 1
            except KeyError:
                has_time = False

            write_var_struct_section(latex, var_struct, struct_name, has_time)

    latex.close()


parser = OptionParser()
parser.add_option("-f", "--file", dest="registry_path",
                  help="Path to Registry file", metavar="FILE")
parser.add_option("-d", "--tex_dir", dest="latex_dir",
                  help="Path to directory with latex addition files.",
                  metavar="DIR")
parser.add_option("-p", "--tex_path", dest="latex_path",
                  help="Path to latex input files that will be written to "
                       "generated latex.", metavar="PATH")

options, args = parser.parse_args()

if not options.registry_path:
    parser.error("Registry file is required")

if not options.latex_dir:
    parser.error('Directory with group latex files is missing. Skipping '
                 'addition of latex files.')
if not options.latex_path:
    parser.error('Need latex path with latex directory.')

latex_missing_string = '{\\bf \color{red} MISSING}'
dimension_table_header = '{\\bf Name} & {\\bf Units} & {\\bf Description}'
variable_table_header = '{\\bf Name} & {\\bf Description}'
namelist_table_header = '{\\bf Name} & {\\bf Description}'

registry_path = options.registry_path

registry_tree = ET.parse(registry_path)

registry = registry_tree.getroot()

write_namelist_input_generated()

# Write file that defines version string for model.
latex = open('define_version.tex', 'w')
try:
    version_string = registry.attrib['version']
except KeyError:
    version_string = '{\\bf MISSING}'
latex.write('\\newcommand{\\version}{%s}\n' % version_string)
latex.close()

write_dimension_table_documentation()

write_namelist_table_documentation()

write_namelist_section_documentation()

write_variable_table_documentation()

write_variable_section_documentation()
