#!/usr/bin/env python
"""
Make a graph file from the MPAS mesh for use in the Metis graph partitioning
software
"""

from mpas_tools.parallel import make_graph_file_entry_point


if __name__ == '__main__':
    make_graph_file_entry_point()

