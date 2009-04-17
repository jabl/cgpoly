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

"""Module for dealing with configuration parameters"""


class CGConfig(object):
    """Base class for objects with config params"""
    
    def __init__(self, conf):
        """Initialize by getting params from file or dict

        Arguments:
        conf: A string containing file name or a dict with configuration values

        """
        if isinstance(conf, str):
            self.c = {}
            execfile(conf, {}, self.c)
        elif not conf:
            self.c = {}
            execfile('cgparams.py', {}, self.c)
        else:
            self.c = conf

    def _get_sigmas(self):
        return (self.c['p_sigma'], self.c['c_sigma'], self.c['i_sigma'])
    sigmas = property(_get_sigmas, None, None, 'Sigmas of the beads')

    def _get_kbt(self):
        return self.c['kb'] * self.c['temp']
    kbt = property(_get_kbt, None, None, 'kb * temp')

    def _get_table_extension(self):
        return self.c['table_extension']
    table_extension = property(_get_table_extension, None, None, 'Table extension')

    def _get_cp_r(self):
        return self.c['cp_r']
    cp_r = property(_get_cp_r, doc='C-P bond length')

    def _get_ip_r(self):
        return self.c['ip_r']
    ip_r = property(_get_ip_r, doc='I-P bond length')

    def _get_avg_r(self):
        return (self.c['ip_r'] + self.c['cp_r']) / 2
    avg_r = property(_get_avg_r, doc='Average bond length')

    def _get_cp_k(self):
        return self.c['cp_k']
    cp_k = property(_get_cp_k, doc='C-P spring constant')

    def _get_ip_k(self):
        return self.c['ip_k']
    ip_k = property(_get_ip_k, doc='I-P spring constant')

    def _get_ipc_theta(self):
        return self.c['ipc_theta']
    ipc_theta = property(_get_ipc_theta, doc='I-P-C equilibrium angle')

    def _get_ipc_k(self):
        return self.c['ipc_k']
    ipc_k = property(_get_ipc_k, doc='I-P-C spring constant')
