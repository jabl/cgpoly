# -*- coding: latin-1 -*-
# Parameters of the coarse-grained BPA-PC model

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
tau_t = 1/(0.5 * time_scaling)


# Harmonic bonds, k and b, see section 4.2.1 in gromacs manual
# Harminic spring constant is taken from gaussian distr
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
c_sigma = 0.349 / sigma
p_sigma = 0.519 / sigma
i_sigma = 0.467 / sigma

# Some technical stuff

# How much tables extend beyond cutoff. This value must match
# table-extension in the .mdp file.
table_extension = 1
