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

"""This module contains functions to generate polymer coordinates."""

from random import uniform
import numpy as np

class ConstraintError(Exception):
    """Exception class for when a coordinate is invalid for some reason."""
    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return repr(self.value)

class Plane(object):
    """Class for representing a plane in 3D space.

    This class is copied from the vasputil.geometry module.

    """

    def __init__(self, points, normal=None):
        """Initialize plane.

        Arguments:
        points: points in the plane.
        normal: Normal vector of the plane.

        If normal is not provided, points must be a sequence or Numpy ndarray
        providing coordinates of 3 points in the plane. If coordinates are
        provided as a Numpy ndarray, each coordinate must be a row vector. If
        normal is provided, a single point is sufficient.

        """
        if normal == None:
            if isinstance(points, np.ndarray) and points.shape != (3,3):
                raise TypeError("Shape of points array must be (3,3).")
            elif len(points) != 3:
                raise TypeError("Points sequence must have 3 elemnents.")
            v1 = points[1] - points[0]
            v2 = points[2] - points[1]
            self.normal = np.cross(v1, v2)
            self.normal /= np.linalg.norm(self.normal)
            self.d_origo = np.dot(self.normal, points[1])
        else:
            self.normal = normal / np.linalg.norm(normal)
            self.d_origo = np.dot(self.normal, points)

    def distance(self, point):
        """Measure the distance between the plane and a point."""
        return abs(np.dot(self.normal, point) - self.d_origo)

# End of class Plane


class WallConstraint(object):
    """A wall constraint"""
    def __init__(self, z=0):
        """Initialize wall constraint

        z: Place wall z coordinate here
        d: Don't let polymer random walk get closer than this
        """
        self.plane = Plane((0,0,z), (0, 0, 1))

    def check(self, coords):
        """Check whether coords are violating constraint
        
        coords: 2x3 array with coordinates in row vectors
        """
        d0 = np.dot(self.plane.normal, coords[0]) - self.plane.d_origo
        d1 = np.dot(self.plane.normal, coords[1]) - self.plane.d_origo
        if np.sign(d0) != np.sign(d1):
            raise ConstraintError('Polymer hit wall')


def fold_coords(coords, box, pbc):
    """Fold coordinates into box or throw exception"""
    for ii, bdim in enumerate(box):
        if not pbc[ii]:
            if coords[ii] < 0. or coords[ii] > bdim:
                raise ConstraintError("Coordinate outside box!")
            continue
        coords[ii] = coords[ii] % bdim
    return coords

def vec_pbc(vec, box, pbc):
    """Vector taking into account periodic boundary conditions.

    This code is adapted from vasputil.geometry module.

    """
    arr = np.array(vec)
    v1 = arr.reshape(arr.size)
    for ii, bdim in enumerate(box):
        if not pbc[ii]:
            continue
        b2 = bdim / 2
        if v1[ii] < -b2:
            v1[ii] += bdim
        elif v1[ii] > b2:
            v1[ii] -= bdim
    return arr


def euler_rot(dr, cc):
    """Create Euler rotation matrix and rotate vector

    This transforms the rotated coordinate system where z is
    aligned on the bead axis back to the usual one.

    dr -- vector that defines the z axis in the rotated coord system
    cc -- Vector in the rotated coord. system to transform

    Returns the vector cc rotated to the usual xyz coord system.

    See http://mathworld.wolfram.com/EulerAngles.html

    We know that in the rotated coordinates the coords of the n-1
    bead is (0,0,0), whereas in the system cartesian coordinates
    they are given by the dr input argument considering the n-1 has
    first been transformed to origo. Hence the rotation matrix can be
    calculated.

    Need to create orthogonal vectors in order to give the rotation
    matrix capability for arbitrary rotations. Or not necessarily
    orthogonal, but at least a set of vectors that span the space.

    """
    dro1 = np.cross(dr, np.array((0, 0, 1)))
    dro2 = np.cross(dr, dro1)
    return np.linalg.solve(np.array((dro1, dro2, dr)), cc)

def gen_chain_rw(nbeads, box, pbc, dist, angle, maxtry=100, constr=[]):
    """Generate a single linear chain using a randow walk

    nbeads -- Number of beads in the chain to generate
    
    box -- Sequence of length 3 containing box dimensions in (x,y,z)
    
    pbc -- Sequence of length 3 describing whether periodic boundaries
           are used.

    dist -- Distance between two beads

    angle -- Angle between two beads in degrees. The angle is defined
             as in Gromacs manual fig 4.7

    maxtry -- Maximum number of tries before giving up

    Returns a numpy array with bead coordinates

    """
    angle = angle * np.pi/180
    coords = np.empty((nbeads, 3))
    # Radius of the cylinder for angle, i.e. how much is the next bead
    # offset from the axis of the two previous beads.
    cyl_rad = dist * np.sin(np.pi - angle)
    # Z coord of current bead where cylinder origo is at the n-1 bead
    cyl_z = dist * np.cos(np.pi - angle)
    for ntry in range(maxtry):
        # Starting coordinates of 1st bead
        for ii, bdim in enumerate(box):
            coords[0, ii] = uniform(0, bdim)
        # Second bead, using spherical coordinates
        theta = uniform(0, 2*np.pi)
        phi = uniform(0, np.pi)
        sinphi = np.sin(phi)
        coords[1] = coords[0] + dist * np.array((
            np.cos(theta) * sinphi,
            np.sin(theta) * sinphi,
            np.cos(phi)))
        try:
            for con in constr:
                con.check(coords[0:2])
            #coords[1] = fold_coords(coords[1], box, pbc)
            # Rest of the chain
            # Use cylindrical coordinates with z axis along n-1 and n-2 beads
            # Origo at the n-1 bead
            for ii in range(2, nbeads):
                theta = uniform(0, 2*np.pi)
                cc = np.array((cyl_rad * np.cos(theta),
                               cyl_rad * np.sin(theta),
                               cyl_z))
                # dr argument to euler_rot must be length 1
                dr = (coords[ii-1] - coords[ii-2]) / dist
                cc = euler_rot(dr, cc)
                # Is this a kludge bugfix, or why doesn't the rotation
                # matrix preserve the length otherwise?
                cc = cc / np.linalg.norm(cc) * dist
                coords[ii] = cc + coords[ii-1]
                # Add checking for constraint errors, going out of the
                # box if not pbc=xyz etc. here.
                for con in constr:
                    con.check(coords[(ii-1):(ii+1)])
        except ConstraintError:
            if ntry == maxtry -1:
                raise ConstraintError("Polymer generation failed, maybe \
                try increasing maxtry or box size?")
            continue # Restart polymer generation
        break
    return coords

def gen_index(mbeads, bpc, nchains):
    """Generate index files for angles

    First generate angles from monomer descriptor. E.g.
    [P,C,P,I] => [PCP, CPI, PIP, IPC]
    Remove duplicates => [PCP, PIP, IPC]
    Generate files => [PCP-angles.ndx, PIP-angles.ndx, IPC-angles.ndx]

    """
    angles = []
    mbl = len(mbeads)
    for ii, bead in enumerate(mbeads):
        jj = (ii + 1) % mbl
        kk = (ii + 2) % mbl
        angles.append([bead + mbeads[jj] + mbeads[kk], ii])

    # Remove duplicates by checking if reversing them matches another
    # angle (insert fancy symmetry group name here to make the
    # crystallography weenies happy)
    ra = []
    for ii, ang in enumerate(angles):
        as = ang[0]
        angs = [a[0] for a in angles[(ii+1):]]
        for jj, ang2 in enumerate(angs):
            if as[::-1] == ang2:
                ra.append((ang, jj + ii))
    for a in ra:
        angles[a[1] + 1].append(a[1]-1)
        angles.remove(a[0])

    # Generate files
    for ang in angles:
        f = open(ang[0] + '-angles.ndx', 'w')
        f.write('[ %s ]\n' % ang[0])
        aind = ang[1:]
        aind.sort()
        for chain in range(nchains):
            for bead in range(bpc-2):
                for a in aind:
                    if bead % mbl == a:
                        bn = bpc*chain + bead 
                        f.write('%i %i %i\n' % (bn+1, bn+2, bn+3))


def gen_bpapc(filename, bpc, nchains, box, pbc, dist, angle, constr=[]):
    """Write out coordinate file for BPA-PC system.

    The file is in GRO format. See Gromacs manual for description.

    """
    f = open(filename, 'w')
    f.write('BPA-PC system with %i chains and %i beads per chain, t = 0.0\n'
            % (nchains, bpc))
    f.write('%i\n' % (bpc * nchains))
    for chain in range(nchains):
        coords = gen_chain_rw(bpc, box, pbc, dist, angle, constr=constr)
        for ii, bead in enumerate(coords):
            if ii % 4 == 0 or ii % 4 == 2:
                btype = 'P'
            elif ii % 4 == 1:
                btype = 'C'
            else:
                btype = 'I'
            btype += str(ii+1)
            f.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"
                    % (chain + 1, 'BPAPC', btype, chain*bpc + ii + 1,
                       bead[0], bead[1], bead[2], 0.0, 0.0, 0.0))
    f.write('%8.4f %8.4f %8.4f\n' % box)
    f.close()
    # Finally generate index files for g_angle
    gen_index(['P', 'C', 'P', 'I'], bpc, nchains)

if __name__ == '__main__':
    from optparse import OptionParser

    usage = """%prog [options] outfile

Generate BPA-PC polymer coordinates, storing them into the outfile in
GRO format. Also generates index files for the angles PCP-angles.ndx,
PIP-angles.ndx, IPC-angles.ndx that can be used with g_angle."""
    
    parser = OptionParser(usage)
    parser.add_option('-b', '--beads-per-chain', dest='bpc', help='Number of beads per chain')
    parser.add_option('-n', '--num-chains', dest='nchains', help='Number of chains')
    parser.add_option('-s', '--size', dest='size', help='Size of simulation cube side')
    parser.add_option('-a', '--angle', dest='angle', help='Angle between beads')
    (options, args) = parser.parse_args()
    if options.bpc:
        bpc = int(options.bpc)
    else:
        bpc = 43
    if options.nchains:
        nchains = int(options.nchains)
    else:
        nchains = 80
    if options.size:
        s = float(options.size)
        box = (s, s, s)
    else:
        box = (7, 7, 7)
    dist = (cp_r + ip_r) / 2
    # Average angle
    # angles in monomer: PCP, IPC, PIP, IPC
    if options.angle:
        angle = float(options.angle)
    else:
        angle = 145.61693134825001
    gen_bpapc(args[0], bpc, nchains, box, (True, True, True), dist, angle)
