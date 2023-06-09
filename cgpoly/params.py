#!/usr/bin/python

# Copyright (C) 2009, 2023 Janne Blomqvist
# SPDX-License-Identifier: Apache-2.0

"""Module for dealing with configuration parameters"""

import os

def find_datadir():
    """Try to find the data files."""
    # Check current directory
    if os.path.exists('cgparams.py'):
        return '.'
    # 1st try didn't succeed. Check if we find the installed data dir
    p = __file__
    while 1:
        p = os.path.dirname(p)
        d = os.path.join(p, 'share/cgpoly')
        if os.path.exists(d):
            return d
        if p == '/':
            return None

def copy_data_files():
    """Copy needed data files to current directory."""
    d = find_datadir()
    if d == '.' or d == None:
        return
    os.system('cp ' + d + '/* .')


class CGConfig(object):
    """Base class for objects with config params"""
    
    def __init__(self, conf=None):
        """Initialize by getting params from file or dict

        Arguments:
        conf: A string containing file name or a dict with configuration values

        """
        if isinstance(conf, str):
            self.c = {}
            execfile(conf, {}, self.c)
        elif not conf:
            self.c = {}
            execfile(os.path.join(find_datadir(), 'cgparams.py'), {}, self.c)
        else:
            self.c = conf

    def _get_sigmas(self):
        return (self.c['p_sigma'], self.c['c_sigma'], self.c['i_sigma'])
    sigmas = property(_get_sigmas, None, None, 'Sigmas of the beads')

    def _get_sigmas_wall(self):
        return (self.c['pw_sigma'], self.c['cw_sigma'], self.c['iw_sigma'])
    sigmas_wall = property(_get_sigmas_wall, None, None, 'Wall sigmas')

    def _get_eps_wall(self):
        return (self.c['pw_eps'], self.c['cw_eps'], self.c['iw_eps'])
    eps_wall = property(_get_eps_wall, None, None, 'Wall epsilons')

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



            
        
