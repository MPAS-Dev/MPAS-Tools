#!/usr/bin/env python
'''
This script marks all of the boundary cells in a domain as Dirichlet velocity
boundaries.
'''

from mpas_tools.landice.boundary import mark_domain_boundaries_dirichlet


if __name__ == '__main__':
    mark_domain_boundaries_dirichlet()

