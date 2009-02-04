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

avg_angle = (exp_pcp + exp_pip + 2*180.)/4

# Table files with non-bonded interactions (WCA potentials)
gen_nb_tables.gen_nb_files()

# Chain topology and bonded interactions
chain_top.genchain_top(bpc, 'bpapc43.itp')

# Finally, coordinates
gen_poly.gen_bpapc('bpapc.gro', bpc, 80, (7,7,7), (True, True, True), dist,
                   avg_angle)

print 'Now create the index file. See readme.txt for details.'
print """Now run grompp and mdrun. For equilibration, there is
equil.mdp, for real runs md.mdp. Equilibration runs use force capping
(soft cores) for nonbonded interactions. E.g.:"""
print 'grompp -f equil.mdp -c bpapc.gro -p bpapc80.top -n index.ndx'
print 'And then mdrun -v'
