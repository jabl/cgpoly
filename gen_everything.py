#!/usr/bin/python

# Copyright (C) 2009 Janne Blomqvist

# Generate everything needed for a coarse-grained BPA-PC system
# with 80 chains and 43 beads per chain in a 7x7x7 nm box

import espresso2gro_table as e2g
import gen_nb_tables
import chain_top
import gen_poly
from params import *

# Beads per chain
bpc = 43

# Average bead-bead distance
dist = (cp_r + ip_r) / 2

# First, convert Espresso format tables for bonded angle potentials to
# Gromacs format
exp_pcp = e2g.espresso2gro_convert_table('PCP_TABLE.DAT', 'table_a0.xvg')
exp_pip = e2g.espresso2gro_convert_table('PIP_TABLE.DAT', 'table_a1.xvg')

# Average IPC angle is 175.7 deg taken from an equilibrated simulation
# (obviously can't never be 180)
avg_angle = (exp_pcp + exp_pip + 2*175.731)/4
# However, this angle produces a too large radius of gyration, so just
# use one that has been determined to produce the correct one
avg_angle = 141

# Table files with non-bonded interactions (WCA potentials)
gen_nb_tables.gen_nb_files()

# Chain topology and bonded interactions
chain_top.genchain_top(bpc, 'bpapc43.itp')

# Finally, coordinates
gen_poly.gen_bpapc('bpapc.gro', bpc, 80, (7,7,7), (True, True, True), dist,
                   avg_angle)
