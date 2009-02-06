#!/usr/bin/python

# Copyright (C) Janne Blomqvist 2009

# Generate tables with the nonbonded interactions for the
# coarse-grained BPA-PC model

# A complication is that we need different sigma_ij for different bead
# types. This can be done with the energygrp_table option in the .mdp
# file. See the thread at
# http://www.mail-archive.com/gmx-users@gromacs.org/msg07992.html

from params import *
import ljpot
import numpy as np
import os

def cutoff_sigma(sig):
    "WCA cutoff from sigma"
    return sig * 2**(1./6)

def maxcut():
    "Max cutoff needed for tables"
    cuts = [cutoff_sigma(x) for x in [c_sigma, p_sigma, i_sigma]]
    # Round up to nearest 0.1 nm, add another 0.1 for safety.
    return np.ceil(max(cuts) * 10) / 10. + table_extension + 0.1

def nb_tables(sigma1, sigma2=None, forcecap=None, dr=0.002):

    """Create tables with pot and force

    The table goes from 0 to cutoff in dr intervals. cutoff is
    determined by sigma. The nonbonded interactions in reduced units
    are according to eq. 5.1 in the Espresso manual with params

    eps = 1
    sigma = reduced per-atomtype sigma
    cutoff = sigma * 2**(1./6) (should be sigma**(1./6) for WCA?)
    shift = 1./4
    offset = 0

    """
    if sigma2 == None:
        sigma2 = sigma1
    sig = (sigma1 + sigma2) / 2
    eps = kb * temp
    cutoff = cutoff_sigma(sig)
    mc = maxcut()
    rr = np.linspace(0, mc, mc/dr)
    pot = np.empty(rr.shape)
    force = pot.copy()
    pot[1:], force[1:] = ljpot.lj(rr[1:], eps, sig, 1./4, cutoff, forcecap)
    pot[0] = pot[1]
    force[0] = force[1]
    return rr, pot, force

def gen_nb_files(forcecap=None):
    """Generate Gromacs nb tables

    The actual data will be as g and g', that is columns 4 and 5. In
    the topology the C6 factor, which should be 1.0, will be
    multiplied with the value taken from the table.

    """
    beads = ['P', 'C', 'I']
    sigmas = [p_sigma, c_sigma, i_sigma]
    for ii, b1 in enumerate(beads):
        for jj, b2 in enumerate(beads[ii:]):
            rr, pot, force = nb_tables(sigmas[ii], sigmas[jj], forcecap)
            fn = 'table_' + b1 + '_' + b2 + '.xvg'
            fout = open(fn, 'w')
            for kk, dist in enumerate(rr):
                fout.write('%5.4f 0 0 %5g %5g 0 0\n'
                           % (dist, pot[kk], force[kk]))
            fout.close()
    # Even though pairwise tables are used for all NB interactions, gromacs
    # complains if it can't find the generic table.xvg used for all others.
    # So just create a dummy.
    os.remove('table.xvg')
    os.symlink('table_P_P.xvg', 'table.xvg')


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--forcecap', dest='fcap', help='Cap maximum force at this value')
    (options, args) = parser.parse_args()
    if options.fcap:
        fcap = float(options.fcap)
    else:
        fcap = None
    gen_nb_files(fcap)
