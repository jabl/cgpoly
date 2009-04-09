#!/usr/bin/python
# -*- coding: latin-1 -*-

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


"""This module contains the stuff to generate a GROMACS topology file
for a single chain, as well as to convert an Espresso angle potential
file to GROMACS format."""

import time, math
import numpy as np
import cgpoly.params

# BPA-PC monomer is [P-C-P-I]

class TopologyGenerator(cgpoly.params.CGConfig):
    """Topology Generator"""

    def genchain_top(self, bpc, fout):
        "Generate topology for BPA-PC chain with bpc beads per chain"
        if isinstance(fout, str):
            fout = open(fout, 'w')
        fout.write("""; BPA-PC chain
; DO NOT EDIT! - This file is automatically generated
; Generated by chain_top.py on %s

[ moleculetype ]
; molname   nrexcl
BPAPC       3

[ atoms ]
; id  type  resnr  residu  atom  cgnr  charge\n""" % time.asctime())
        for b in range(1, bpc+1):
            bn = (b-1) % 4
            bead = ['P', 'C', 'P', 'I'][bn]
            # For bead name just use bead type + monomer number
            bname = bead + str(b)
            fout.write("%4i     %s      1   BPAPC   %3s  %4i     0.0\n"
                       % (b, bead, bname, b))
        fout.write("""
[ bonds ]
; i j func length force\n""")
        for b in range(1, bpc):
            bn = (b-1) % 4
            blen = [self.cp_r, self.cp_r, self.ip_r, self.ip_r][bn]
            kforce = [self.cp_k, self.cp_k, self.ip_k, self.ip_k][bn]
            fout.write("%4i  %4i  1  %7f  %7f\n" % (b, b+1, blen, kforce))

        fout.write("""
[ angles ]
; i j k func angle(deg) force\n""")
        for b in range(1, bpc-1):
            bn = (b-1) % 4
            if bn == 1 or bn == 3:
                ban = self.ipc_theta
                force = self.ipc_k
                fout.write("%4i  %4i  %4i  1  %7.4f %9.4f\n"
                           % (b, b+1, b+2, ban, force))
            elif bn == 0 or bn == 2:
                tablenum = bn/2
                fout.write("%4i  %4i  %4i  8  %i 1.0\n" % (b, b+1, b+2, tablenum))

        
    def espresso2gro_convert_table(self, fin, fout, plot=False):
        """Convert angle potential table from Espresso format to Gromacs"""
        if isinstance(fin, str):
            fin = open(fin)
        if isinstance(fout, str):
            fout = open(fout, 'w')
        pots = []
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
            else:
                angle, force, pot = [float(x) for x in line.split()]
                angle = angle * 180. / math.pi
                pot *= self.kbt
                pots.append((angle, pot, -force))

        # Espresso tables don't contain the 0.0 angle, insert it manually
        # It's practically never going to happen, so the actual values
        # don't matter that much. I.e. don't bother with some fancy
        # extrapolation scheme.
        pots.insert(0, (0.0, pots[0][1], pots[0][2]))
    
        pots = np.asarray(pots)
        for line in pots:
            # As the forces in Espresso tables are suspiciously noisy,
            # just let Gromacs use forces generated by taking the
            # numerical derivative of the potential (set force=0 in table)
            fout.write("%20.10f  %20.12f  %20.12f\n" % (line[0], line[1], 0.0))
    
        # Calculate expectation value from distribution
        # Distribution is Boltzmann inversion of potential
    
        # Ignore the partition function Q, we're just interested in
        # the shape of the distribution
    
        w = np.exp(-pots[:,1] / self.kbt)
        norm = np.sum(w)
        kernel = pots[:,0] * w / norm
        eval = np.sum(kernel)
        if plot:
            import pylab
            pylab.plot(pots[:,0], pots[:,1], label='Potential')
            pylab.plot(pots[:,0], pots[:,2], label='Force')
            pylab.plot(pots[:,0], w*50, label='Distribution (scaled x50)')
            pylab.legend()
            pylab.xlim(0, 180)
            pylab.show()
        return eval


if __name__ == '__main__':
    from optparse import OptionParser
    usage = """%prog [options] outfile

Generate gromacs topology file with bonded interaction
parameters. Alternatively, if the -e option is given, convert an
expresso format angle potential table to GROMACS format."""
    parser = OptionParser(usage)
    parser.add_option('-b', '--beads-per-chain', dest='bpc', help='Generate topology for chain with BPC number of beads')
    parser.add_option('-e', '--espresso-table', dest='etab', help='Convert Espresso table to GROMACS format')
    (options, args) = parser.parse_args()
    if options.bpc:
        bpc = int(options.bpc)
    else:
        bpc = 43
    if options.etab:
        eval = espresso2gro_convert_table(options.etab, args[0], True)
        print 'Expectation value: ', eval
    else:
        genchain_top(bpc, args[0])
