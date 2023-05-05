#!/usr/bin/python

# Copyright (C) 2009, 2023 Janne Blomqvist
# SPDX-License-Identifier: Apache-2.0

"""Generate tables with the nonbonded interactions for the
coarse-grained BPA-PC model

A complication is that we need different sigma_ij for different bead
types. This can be done with the energygrp_table option in the .mdp
file. See the thread at
http://www.mail-archive.com/gmx-users@gromacs.org/msg07992.html"""

import ljpot
import numpy as np
import os
import cgpoly.params as cp

def cutoff_sigma(sig):
    "WCA cutoff from sigma"
    return sig * 2**(1./6)

class NonBonded(cp.CGConfig):
    """Functionality for creating tables for nonbonded interactions"""

    def __init__(self, conf=None, wall=False):
        """Initialize nonbonded table creator"""
        super(NonBonded, self).__init__(conf)
        self.wall = wall

    def cutoff_wall(self, sig):
        """Cutoff for wall 10-4 LJ pot"""
        return sig * self.c['cutoff_wall']
    
    def maxcut(self):
        "Max cutoff needed for tables"
        cuts = [cutoff_sigma(x) for x in self.sigmas]
        if self.wall:
            for ii, x in enumerate(self.sigmas_wall):
                if self.eps_wall[ii] != None:
                    cuts.append(self.cutoff_wall(x))
                else:
                    cuts.append(cutoff_sigma(x))
        # Round up to nearest 0.1 nm, add another 0.1 for safety.
        return np.ceil(max(cuts) * 10) / 10. + self.table_extension + 0.1

    def nb_tables(self, sigma1, sigma2=None, forcecap=None, dr=0.002,
                  wall=False, eps=None):

        """Create tables with pot and force

        The table goes from 0 to cutoff in dr intervals. cutoff is
        determined by sigma. The nonbonded interactions in reduced units
        are according to eq. 5.1 in the Espresso manual with params

        eps = 1
        sigma = reduced per-atomtype sigma
        cutoff = sigma * 2**(1./6) (should be sigma**(1./6) for WCA?)
        shift = 1./4
        offset = 0

        if wall=True, then generates bead-wall 10-4 LJ tables instead.
        
        """
        if sigma2 == None:
            sigma2 = sigma1
        sig = (sigma1 + sigma2) / 2
        if wall and eps != None:
            cutoff = self.cutoff_wall(sig)
        else:
            cutoff = cutoff_sigma(sig)
        mc = self.maxcut()
        rr = np.arange(0, mc + dr, dr)
        pot = np.empty(rr.shape)
        force = pot.copy()
        if wall and eps != None:
            #pot[1:], force[1:] = ljpot.ljwall(rr[1:], eps, sig, 3./5, cutoff, forcecap)
            pot[1:], force[1:] = ljpot.ljwall(rr[1:], eps, sig, 0., None, forcecap)
        else:
            if eps == None:
                eps = self.kbt
            pot[1:], force[1:] = ljpot.lj(rr[1:], eps, sig, 1./4, cutoff, forcecap)
        pot[0] = pot[1]
        force[0] = force[1]
        return rr, pot, force

    def gen_nb_files(self, forcecap=None):
        """Generate Gromacs nb tables
        
        The actual data will be as g and g', that is columns 4 and 5. In
        the topology the C6 factor, which should be 1.0, will be
        multiplied with the value taken from the table.
        
        """
        beads = ['P', 'C', 'I']
        sigmas = self.sigmas
        for ii, b1 in enumerate(beads):
            for jj, b2 in enumerate(beads[ii:]):
                rr, pot, force = self.nb_tables(sigmas[ii], sigmas[jj + ii],
                                                forcecap)
                fn = 'table_' + b1 + '_' + b2 + '.xvg'
                fout = open(fn, 'w')
                for kk, dist in enumerate(rr):
                    fout.write('%5.4f 0 0 %5g %5g 0 0\n'
                               % (dist, pot[kk], force[kk]))
                fout.close()
        # If wall potentials are also desired, generate them
        if self.wall:
            for ii, bb in enumerate(beads):
                rr, pot, force = self.nb_tables(self.sigmas_wall[ii],
                                                forcecap=forcecap, wall=True,
                                                eps=self.eps_wall[ii])
                fnb = 'table_' + bb
                fn0 = fnb + '_wall0.xvg'
                fout = open(fn0, 'w')
                for kk, dist in enumerate(rr):
                    fout.write('%5.4f 0 0 %5g %5g 0 0\n'
                               % (dist, pot[kk], force[kk]))
                fout.close()
                fn1 = fnb + '_wall1.xvg'
                if os.path.exists(fn1):
                    os.remove(fn1)
                os.symlink(fn0, fn1)
        # Even though pairwise tables are used for all NB interactions, gromacs
        # complains if it can't find the generic table.xvg used for all others.
        # So just create a dummy.
        if os.path.exists('table.xvg'):
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
    nb = NonBonded()
    nb.gen_nb_files(fcap)
