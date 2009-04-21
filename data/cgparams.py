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


"""Parameters of the coarse-grained BPA-PC model

See
Abrams, Macromolecules 36, 260-267 (2003)"""

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

# System scaling factor for reduced units
rsigma = 4.41
# For normal units, for reduced units set sigma = rsigma
sigma = 1.

# Timestep, in Espresso with reduced units 0.001
time_scaling = rsigma / math.sqrt(kb*temp / 73)
timestep = 0.001 * time_scaling
# gamma = 1/tau
tau_t = 1/0.5 * time_scaling


# Harmonic bonds, k and b, see section 4.2.1 in gromacs manual
# Harmonic spring constant is taken from gaussian distr
# 1/(stddev/sigma)**2
cp_k = (sigma / 0.04) ** 2
cp_r = 0.356 / sigma

ip_k = cp_k
ip_r = 0.293 / sigma

# Angle for IPC, very stiff 180 deg
ipc_theta = 180.0
# Gromacs wants angle force in kJ/mol/rad**2
# The force constant is 0.1 kBT/deg**2
ipc_k = 0.1 * kb * temp * (180. / math.pi)**2

# Non-bonded parameters
# Note the value of c_sigma in Abrams MM 36 260 is incorrect, this is the
# correct value!
c_sigma = 0.349 / sigma
p_sigma = 0.467 / sigma
i_sigma = 0.519 / sigma


# Wall interaction 10-4 LJ parameters
ev2kjmol = 96.48538 # 1 eV = 96.48538 kJ/mol
cw_eps = 1.89 * ev2kjmol
pw_eps = 0.66 * ev2kjmol
iw_eps = 0.
cw_sigma = 2.97948739942 / 10 / sigma
pw_sigma = 3.42799564404 / 10 / sigma
iw_sigma = 4. / 10 / sigma # Purely repulsive

# Some technical stuff

# How much tables extend beyond cutoff. This value must match
# table-extension in the .mdp file.
table_extension = 1
