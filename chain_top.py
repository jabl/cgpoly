#!/usr/bin/python
# -*- coding: latin-1 -*-

# Copyright Janne Blomqvist (C) 2009

# This file contains the parameters for the bonded interactions in the
# coarse-grained BPA-PC model, and is used to generate a GROMACS
# topology file for a single chain.

import time, math
from params import *

# BPA-PC monomer is [P-C-P-I]


def genchain_top(bpc, fout):
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
        blen = [cp_r, cp_r, ip_r, ip_r][bn]
        kforce = [cp_k, cp_k, ip_k, ip_k][bn]
        fout.write("%4i  %4i  1  %7f  %7f\n" % (b, b+1, blen, kforce))

    fout.write("""
[ angles ]
; i j k func angle(deg) force\n""")
    for b in range(1, bpc-1):
        bn = (b-1) % 4
        if bn == 1 or bn == 3:
            ban = ipc_theta
            force = ipc_k
            fout.write("%4i  %4i  %4i  1  %7.4f %9.4f\n"
                       % (b, b+1, b+2, ban, force))
        elif bn == 0 or bn == 2:
            tablenum = bn/2
            fout.write("%4i  %4i  %4i  8  %i 1.0\n" % (b, b+1, b+2, tablenum))

if __name__ == '__main__':
    from optparse import OptionParser
    usage = """%prog [options] outfile"""
    parser = OptionParser(usage)
    parser.add_option('-b', '--beads-per-chain', dest='bpc', help='Generate topology for chain with BPC number of beads')
    (options, args) = parser.parse_args()
    if options.bpc:
        bpc = int(options.bpc)
    else:
        bpc = 43
    genchain_top(bpc, args[0])
