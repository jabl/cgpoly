#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2009, 2010, 2023 Janne Blomqvist
# SPDX-License-Identifier: Apache-2.0

"""Different variants of LJ potentials"""

import numpy as np

def ljwall_scalar(r, epsilon, sigma, cutoff, shift, offset):
    """Helper routine assuming scalar arguments."""
    sf = sigma / (r - offset)
    if r < cutoff:
        return 2 * pi * epsilon * sigma ** 2 \
               * ((2./5) * sf**10 - sf**4 + shift)
    else:
        return 0.

def ljwallold(r, epsilon=1., sigma=1., cutoff=1., shift=3./5, offset=0.):
    """10-4 LJ potential."""

    if isinstance(r, np.ndarray):
        pot = empty(r.size)
        for ind, rr in enumerate(r):
            pot[ind] = ljwall_scalar(rr, epsilon, sigma, cutoff, shift, offset)
        return pot
    else:
        return ljwall_scalar(r, epsilon, sigma, cutoff, shift, offset)

def ljwall(r, epsilon=1., sigma=1., shift=3./5, cutoff=None, forcecap=None):
    """10-4 LJ potential for walls

    Since in this case parameters are taken from DFT calculations with
    surface interaction, the 2*pi prefactor is already included in the
    epsilon.

    However, the max potential, which is at r=sigma, is not -epsilon
    as in regular 12-6 LJ potential, but rather (-3/5)*epsilon. Hence
    we must take this into account.

    """
    sf4 = (sigma / r)**4
    sf10 = (sigma / r)**10
    #pot = 2 * np.pi * epsilon * ( (2. / 5) * sf10 - sf4 + shift)
    #force = 2 * np.pi * epsilon * ( (2.*10 / 5) * sf10 - 4 * sf4) / r
    pot = (5. / 3) * epsilon * ( (2. / 5) * sf10 - sf4 + shift)
    force = (5. / 3) * epsilon * ( (2.*10 / 5) * sf10 - 4 * sf4) / r
    if cutoff != None and isinstance(pot, np.ndarray):
        pot = np.where(r <= cutoff, pot, 0.)
        force = np.where(r <= cutoff, force, 0.)
    if forcecap != None:
        # Force capped potential according to Auhl et al
        # JCP 119, 12718 (2003)
        # First find the cutoff radius for force capping
        ind = (force > forcecap).sum() - 1
        rcut = r[ind]
        potcut, fcut = lj(rcut, epsilon, sigma, shift, cutoff)
        # Extra - on the first term since lj() return -force
        pot = np.where(force <= forcecap, pot, (rcut - r) * fcut + potcut)
        force = np.where(force <= forcecap, force, forcecap)
    return pot, force


def lj(r, epsilon=1., sigma=1., shift=0, cutoff=None, forcecap=None):
    """Standard 12-6 LJ potential"""
    sf6 = (sigma / r)**6
    pot = 4 * epsilon * (sf6**2 - sf6 + shift)
    force = 4 * epsilon * (12 * sf6**2 - 6 * sf6) / r
    if cutoff != None and isinstance(pot, np.ndarray):
        pot = np.where(r <= cutoff, pot, 0.)
        force = np.where(r <= cutoff, force, 0.)
    if forcecap != None:
        # Force capped potential according to Auhl et al
        # JCP 119, 12718 (2003)
        # First find the cutoff radius for force capping
        ind = (force > forcecap).sum() - 1
        rcut = r[ind]
        potcut, fcut = lj(rcut, epsilon, sigma, shift, cutoff)
        # Extra - on the first term since lj() return -force
        pot = np.where(force <= forcecap, pot, (rcut - r) * fcut + potcut)
        force = np.where(force <= forcecap, force, forcecap)
    return pot, force

if __name__ == '__main__':
    from pylab import *
    import cgpoly.params
    import os.path
    ll = linspace(0.2, 4.0, 1000)
    #plot(ll, ljwall(ll), label='LJ Wall default')
    #plot(ll, ljwall(ll, shift=0.), label='ljwall shift=0')
    plot(ll, lj(ll)[0], label='12-6 LJ sigma=1, epsilon=1')
    plot(ll, lj(ll)[1], label='12-6 LJ force sigma=1, epsilon=1')
    # C-C repulsion, reduced
    execfile(os.path.join(cgpoly.params.find_datadir(), 'cgparams.py'))
    rsigma = 1 # Normal units
    plot(ll, lj(ll, 1, c_sigma/rsigma, 0.25, 2**(1./6)*c_sigma/rsigma)[0], label='12-6 WCA C')
    plot(ll, lj(ll, 1, c_sigma/rsigma, 0.25, 2**(1./6)*c_sigma/rsigma)[1], label='12-6 force WCA C')
    # P-P repulsion, reduced
    plot(ll, lj(ll, 1, p_sigma/rsigma, 0.25, 2**(1./6)*p_sigma/rsigma)[0], label='12-6 WCA P')
    #Carbonate wall
##     c_sigma = 1 # scaled
##     plot(ll, ljwall(ll, cutoff=3*c_sigma, shift=0.), label='Carbonic acid')
    
##     # Phenol
##     p_sigma = 1 # scaled, again
##     plot(ll, ljwall(ll, cutoff=3*p_sigma, shift=0.), label='Phenol')
##     # Isopropylene, purely repulsive
##     i_sigma = 1
##     plot(ll, ljwall(ll, cutoff=i_sigma), label='Isopropylene (pure repulsive)')
    ylim(-4, 2)
    legend()
    show()
