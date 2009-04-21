#!/usr/bin/python

# Copyright (C) 2009 Janne Blomqvist

#  This file is part of cgpoly.

#  cgpoly is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.

#  cgpoly is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with cgpoly.  If not, see <http://www.gnu.org/licenses/>.


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

    def __init__(self, conf):
        """Initialize nonbonded table creator"""
        super(NonBonded, self).__init__(conf)
        
    def maxcut(self):
        "Max cutoff needed for tables"
        cuts = [cutoff_sigma(x) for x in self.sigmas]
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
        cutoff = cutoff_sigma(sig)
        mc = self.maxcut()
        rr = np.linspace(0, mc, mc/dr)
        pot = np.empty(rr.shape)
        force = pot.copy()
        if wall:
            #pot[1:], force[1:] = ljpot.ljwall(rr[1:], eps, sig, 3./5, cutoff, forcecap)
            pot[1:], force[1:] = ljpot.ljwall(rr[1:], eps, sig, 0., cutoff, forcecap)
        else:
            if eps == None:
                eps = self.kbt
            pot[1:], force[1:] = ljpot.lj(rr[1:], eps, sig, 1./4, cutoff, forcecap)
        pot[0] = pot[1]
        force[0] = force[1]
        return rr, pot, force

    def gen_nb_files(self, forcecap=None, wall=False):
        """Generate Gromacs nb tables
        
        The actual data will be as g and g', that is columns 4 and 5. In
        the topology the C6 factor, which should be 1.0, will be
        multiplied with the value taken from the table.
        
        """
        beads = ['P', 'C', 'I']
        sigmas = self.sigmas
        for ii, b1 in enumerate(beads):
            for jj, b2 in enumerate(beads[ii:]):
                rr, pot, force = self.nb_tables(sigmas[ii], sigmas[jj], forcecap)
                fn = 'table_' + b1 + '_' + b2 + '.xvg'
                fout = open(fn, 'w')
                for kk, dist in enumerate(rr):
                    fout.write('%5.4f 0 0 %5g %5g 0 0\n'
                               % (dist, pot[kk], force[kk]))
                fout.close()
        # If wall potentials are also desired, generate them
        if wall:
            for ii, bb in enumerate(beads):
                rr, pot, force = self.nb_tables(self.sigmas_wall[ii],
                                                forcecap=forcecap, wall=wall,
                                                eps=self.eps_wall[ii])
                fn = 'table_' + bb + '_W.xvg'
                fout = open(fn, 'w')
                for kk, dist in enumerate(rr):
                    fout.write('%5.4f 0 0 %5g %5g 0 0\n'
                               % (dist, pot[kk], force[kk]))
                fout.close()
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
    gen_nb_files(fcap)
