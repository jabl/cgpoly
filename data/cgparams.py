# -*- coding: utf-8 -*-

# Copyright (C) 2009, 2011 Janne Blomqvist

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


"""Parameters of the coarse-grained BPA-PC model

See
Abrams, Macromolecules 36, 260-267 (2003)


Wall interaction parameters from 

Blomqvist and Salo, Physical Review B 84, 153410 (2011)
"""

import math

# Remember, Espresso uses Å, Gromacs uses nm!
# Of course, scaled by sigma

# In Gromacs, Boltzmann constant is
# kB = 0.008 314 51 kJ / mol / K
kb = 0.00831451
# Due to reduced units, T = 1
# In normal units
temp = 570 # K

# First, BPA-PC parameters

# System scaling factor for reduced units in Espresso [Å/sigma]
_rsigma = 4.41
# For normal units, for reduced units set sigma = rsigma
#sigma = 1.

# Timestep, in Espresso with reduced units 0.001
#time_scaling = rsigma * math.sqrt(73 / (kb*temp))
# Factor of 10 is due to rsigma being in units of Å=> 1e-10, and then 
# a factor of 1e9 due to GROMACS time unit being ps.
time_scaling = _rsigma * math.sqrt(73) / 10
timestep = 0.001 * time_scaling
# gamma = 1/tau
tau_t = 1/0.5 * time_scaling


# Harmonic bonds, k and b, see section 4.2.1 in gromacs manual
# Harmonic spring constant is taken from gaussian distr
# 1/(stddev/sigma)**2
# std.dev is 0.04 Å = 0.004 nm
cp_k = (1 / 0.004) ** 2
cp_r = 0.356

ip_k = cp_k
ip_r = 0.293

# Angle for IPC, very stiff 180 deg
ipc_theta = 180.0
# Gromacs wants angle force in kJ/mol/rad**2
# The force constant is 0.1 kBT/deg**2
ipc_k = 0.1 * kb * temp * (180. / math.pi)**2

# Non-bonded parameters
# Note the value of c_sigma in Abrams MM 36 260 is incorrect, this is the
# correct value!
c_sigma = 0.349
p_sigma = 0.467
i_sigma = 0.519


# Wall interaction 10-4 LJ parameters
ev2kjmol = 96.48538 # 1 eV = 96.48538 kJ/mol

cw_eps = 0.57 * ev2kjmol
pw_eps = 0.50 * ev2kjmol
iw_eps = 0.50 * ev2kjmol 
cw_sigma = 3.12 / 10 # Factor of 10 due to Angstrom -> nm
pw_sigma = 3.02 / 10
iw_sigma = 3.43 / 10

cutoff_wall = 3 # Multiply sigma with this to get the cutoff range

# Some technical stuff

# How much tables extend beyond cutoff. This value must match
# table-extension in the .mdp file.
table_extension = 1
