#!/usr/bin/python
from optparse import OptionParser
import xml.etree.ElementTree as ET
import os


parser = OptionParser()
parser.add_option("--major", action="store_true", dest="major", help="Increment Major Version (Auto-resets minor version.")
parser.add_option("--minor", action="store_true", dest="minor", help="Increment Minor Version.")

options, args = parser.parse_args()

for r, d, f in os.walk("."):
    for files in f:
        if files.endswith(".xml"):
            path = os.path.join(r, files)
            registry_tree = ET.parse(path)
            registry = registry_tree.getroot()
            version = registry.attrib['version']
            version = version.split('.')
            major_ver = int(version[0])
            minor_ver = int(version[1])

            if options.major:
                major_ver = major_ver + 1
                minor_ver = 0
            elif options.minor:
                minor_ver = minor_ver + 1

            print "%s version: %d.%d"%(path, major_ver, minor_ver)

            registry.set("version", "%d.%d"%(major_ver, minor_ver))
            registry_tree.write(path)
