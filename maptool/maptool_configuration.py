#!/usr/bin/env python
from   __future__ import division, unicode_literals, print_function
from   pymatgen import Structure,Molecule,Lattice,Element
from   itertools import combinations
from   maptool.mathematics import length_vectors, angle_vectors, wrap2_pmhalf,\
       unit_vector, rotation_matrix_around_axis_angle, angle_vector
from   math import ceil, sqrt, cos, sin, radians, acos
import numpy as np
import itertools
import random
import sys

eps=1e-10
two_pi=2*np.pi
class MLattice(Lattice):

    def align_with_axis(self, axis=0, round_decimals=14):
        a = self.matrix[0]
        if axis == 0:
            b = np.array([1, 0, 0])
        elif axis == 1:
            b = np.array([0, 1, 0])
        elif axis == 2:
            b = np.array([0, 0, 1])
        else:
            raise ValueError('Axis must be an integer in (0,1,2)')
        if np.linalg.norm(np.cross(a, b)) < eps:
            return
        c = unit_vector(np.cross(a, b))
        av = angle_vector(a, b)
        rotation_matrix = rotation_matrix_around_axis_angle(c, av)
        return np.dot(rotation_matrix, self.matrix.T).T.round(round_decimals)

    def align_with_plane(self, axis=2, round_decimals=14):
        a = self.matrix[1]

        if axis == 0:
            p1 = np.array([0, 1, 0])
            p2 = np.array([0, 0, 1])
        elif axis == 1:
            p1 = np.array([0, 0, 1])
            p2 = np.array([1, 0, 0])
        elif axis == 2:
            p1 = np.array([1, 0, 0])
            p2 = np.array([0, 1, 0])
        else:
            raise ValueError('Axis must be an integer in (0,1,2)')

        if np.linalg.norm(np.cross(p1, a)) < 1E-10:
            return
        c = unit_vector(np.cross(p1, a))
        # print c
        # A vector perpendicular to the plane
        vector_plane = unit_vector(np.cross(p1, p2))
        v1_u = unit_vector(vector_plane)
        v2_u = unit_vector(c)
        proj = np.dot(v1_u, v2_u)
        # print 'Projection', proj
        av = np.arccos(proj)
        # import math
        # print 'Angle=', math.degrees(av)
        rotation_matrix = rotation_matrix_around_axis_angle(p1, -av)
        cell = np.dot(rotation_matrix, self.matrix.T).T.round(round_decimals)
        # print '-->',cell[1], vector_plane
        # print 'PROJECTION', np.dot(cell[1], vector_plane)
        if np.abs(np.dot(cell[1], vector_plane)) > 1E-10:
            # print 'Failed projection', np.dot(cell[1], vector_plane)
            # print cell
            rotation_matrix = rotation_matrix_around_axis_angle(p1, av)
            cell = np.dot(rotation_matrix, self.matrix.T).T.round(round_decimals)
            if np.dot(cell[1], vector_plane) > 1E-10:
                # print 'Failed projection', np.dot(cell[1], vector_plane)
                # print cell
                pass
        return cell

    def minimal_distances(self, reduced1, reduced2):
        """
        Computes a matrix with the minimal distances between
        two sets of points represented as reciprocal coordinates

        :param reduced1: List or array of reduced coordinates for the first
                            set of points
        :param reduced2: Second set of points

        Example:
        >>> import numpy as np
        >>> lattice = Lattice.random_cell('C2')
        >>> r1 = np.random.rand(4, 3)
        >>> r2 = np.random.rand(4, 3)
        >>> dist, close_imgs = lattice.minimal_distances(r1, r2)

        >>> solution = np.zeros((len(r1), len(r2)))

        >>> for i in range(len(r1)):
        ...                reduced2 = r2[j]+close_imgs[i, j]
        ...                cartesian1 = lattice.reduced2cartesian(reduced1)
        ...                cartesian2 = lattice.reduced2cartesian(reduced2)
        ...                diff_vector = cartesian2 - cartesian1
        ...                solution[i, j] = np.linalg.norm(diff_vector)

        >>> solution - dist < 1E-5
        array([[ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True],
        >>> for i in range(len(r1)):
        ...        for j in range(len(r2)):
        ...                reduced1 = r1[i]
        ...                reduced2 = r2[j]+close_imgs[i, j]
        ...                cartesian1 = lattice.reduced2cartesian(reduced1)
        ...                cartesian2 = lattice.reduced2cartesian(reduced2)
        ...                diff_vector = cartesian2 - cartesian1
        ...                solution[i, j] = np.linalg.norm(diff_vector)

        >>> solution - dist < 1E-5
        array([[ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True]], dtype=bool)

        """
        # Just in case of one single coordinate
        reduced1, reduced2 = np.atleast_2d(reduced1, reduced2)

        images = np.array([list(i) for i in itertools.product([-1, 0, 1], repeat=3)])

        red2_images = reduced2[:, None, :] + images[None, :, :]

        cartesian1 = self.reduced2cartesian(reduced1)
        cartesian2 = self.reduced2cartesian(red2_images)

        diff_vectors = cartesian2[None, :, :, :] - cartesian1[:, None, None, :]

        distances = np.sum(diff_vectors * diff_vectors, axis=3)

        close_images = np.zeros((len(reduced1), len(reduced2), 3))
        for i in range(len(reduced1)):
            for j in range(len(reduced2)):
                dij = distances[i, j]
                close_images[i, j] = images[dij == min(dij)][0]

        return np.min(distances, axis=2) ** 0.5, close_images

    def distance2(self, x1, x2, option='reduced', radius=20, limits=None):

        # Compute the vector from x1 to x2
        dv = np.array(x2) - np.array(x1)

        # If we are not in reduced coordinates,
        # Convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)

        if limits is None:
            limits = np.zeros(3, dtype=int)
            corners = self.get_wigner_seitz_container()

            limits[0] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][0] for x in corners]))))), 5) # 5 compare 2pi 
            limits[1] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][1] for x in corners]))))), 5)
            limits[2] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][2] for x in corners]))))), 5)
 # ceil ---- diagonal !!!!
        ret = {}
        for i0 in np.arange(-limits[0], limits[0] + 1):
            for i1 in np.arange(-limits[1], limits[1] + 1):
                for i2 in np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + np.array([i0, i1, i2])
                    norm2 = np.dot(np.dot(dtot, self.metric_tensor), dtot)
                    if norm2 < radius * radius:
                        ret[(i0, i1, i2)] = {'distance': sqrt(norm2), 'image': dtot}

        return ret

    def get_wigner_seitz_container(self):
        """
        Compute the corners of the box that contains the Wigner-Seitz cell
        Note that this is not the standard reciprocal lattice used for 
        solid state physics with a factor of 2 *pi. 

        :return: dict : dictionary with values numpy arrays 
        """
        ret = {}
        for i in itertools.product((-1, 1), repeat=3):
            rec=self.reciprocal_lattice.matrix/two_pi
            rec_metric=np.dot(rec,rec.T)
            ret[i] = np.dot(rec_metric, i * np.diagonal(self.metric_tensor)) # reduce time !!!
        return ret

    def minimal_distance(self, x1, x2, option='reduced'):
        distances_dict = self.distance2(x1, x2, option=option)
        mindist = sys.float_info.max
        for k in distances_dict:
            if 0 < distances_dict[k]['distance'] < mindist:
                mindist = distances_dict[k]['distance']
        return mindist

    def reduced2cartesian(self, x):
        return np.dot(x, self.matrix)
    
    def cartesian2reduced(self,x):
        return np.dot(x,self.inv_matrix)
           
    def distances_in_sphere(self, x1, x2, radius, option='reduced', exclude_out_sphere=True, sort_by_distance=True):
        """
        Returns all the distances between two positions x1 and x2
        taking into account the periodicity of the cell

        :param sort_by_distance:
        :param exclude_out_sphere:
        :param x1:
        :param x2:
        :param radius:
        :param option:
        :return:
        """
        # Compute the vector from x1 to x2
        dv = np.array(x2) - np.array(x1)

        # If the positions are not in reduced coordinates,convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)
        # log.debug('The wrap vector is: %7.3f %7.3f %7.3f' % tuple(dwrap))

        # We need to compute the distances between equivalent faces
        # For that we to compute the perpendicular distance between
        # two faces, the reciprocal lattice produces the inverse of
        # those distances
        recp_len = np.array(self.reciprocal_lattice.abc)/two_pi
        # The reciprocal (2*pi) is not necessary
        limits = np.ceil(radius * recp_len).astype(int)
        # log.debug('The limits are: %d %d %d' % tuple(limits))

        ret = {}
        # The total size is the product of the number in the 3 directions
        # that is double the limits plus 1 to include the zero
        total_size = np.prod(2 * limits + 1)
        # log.debug('The total size is: %d' % total_size)

        ret['distance'] = np.zeros(total_size)
        ret['image'] = np.zeros((total_size, 3), dtype=int)
        ret['dwrap'] = dwrap
        ret['limits'] = limits

        index = 0
        for i0 in np.arange(-limits[0], limits[0] + 1):
            for i1 in np.arange(-limits[1], limits[1] + 1):
                for i2 in np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + np.array([i0, i1, i2])
                    norm2 = np.dot(np.dot(dtot, self.metric_tensor), dtot)
                    ret['distance'][index] = sqrt(norm2)
                    ret['image'][index] = np.array([i0, i1, i2])
                    index += 1

        # Exclude distances out of sphere
        if exclude_out_sphere:
            inside_sphere = ret['distance'] <= radius
            if sum(inside_sphere) > 0:
                ret['distance'] = ret['distance'][inside_sphere]
                ret['image'] = ret['image'][inside_sphere]

        if sort_by_distance:
            sorted_indices = np.argsort(ret['distance'])
            ret['distance'] = ret['distance'][sorted_indices]
            ret['image'] = ret['image'][sorted_indices]

        return ret

    def minimal_distances(self, reduced1, reduced2):
        """
        Computes a matrix with the minimal distances between
        two sets of points represented as reciprocal coordinates

        :param reduced1: List or array of reduced coordinates for the first
                            set of points
        :param reduced2: Second set of points

        Example:
        >>> import numpy as np
        >>> lattice = Lattice.random_cell('C2')
        >>> r1 = np.random.rand(4, 3)
        >>> r2 = np.random.rand(4, 3)
        >>> dist, close_imgs = lattice.minimal_distances(r1, r2)

        >>> solution = np.zeros((len(r1), len(r2)))

        >>> for i in range(len(r1)):
        ...                reduced2 = r2[j]+close_imgs[i, j]
        ...                cartesian1 = lattice.reduced2cartesian(reduced1)
        ...                cartesian2 = lattice.reduced2cartesian(reduced2)
        ...                diff_vector = cartesian2 - cartesian1
        ...                solution[i, j] = np.linalg.norm(diff_vector)

        >>> solution - dist < 1E-5
        array([[ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True]], dtype=bool)

        """
        # Just in case of one single coordinate
        reduced1, reduced2 = np.atleast_2d(reduced1, reduced2)

        images = np.array([list(i) for i in itertools.product([-1, 0, 1], repeat=3)])

        red2_images = reduced2[:, None, :] + images[None, :, :]

        cartesian1 = self.reduced2cartesian(reduced1)
        cartesian2 = self.reduced2cartesian(red2_images)

        diff_vectors = cartesian2[None, :, :, :] - cartesian1[:, None, None, :]

        distances = np.sum(diff_vectors * diff_vectors, axis=3)

        close_images = np.zeros((len(reduced1), len(reduced2), 3))
        for i in range(len(reduced1)):
            for j in range(len(reduced2)):
                dij = distances[i, j]
                close_images[i, j] = images[dij == min(dij)][0]

        return np.min(distances, axis=2) ** 0.5, close_images

class MStructure(Structure):

    @property
    def is_perfect(self):
       return  int(self.composition.num_atoms)==self.num_sites and min(self.composition.as_dict().values()) >= 1.0 

    @property
    def natom(self):
#        assert self.composition.num_atoms - int(self.composition.num_atoms) < eps
        assert self.is_perfect
        return int(self.composition.num_atoms)

    def sort_sites_using_list(self, sorted_indices):
        sorted_indices = np.array([int(x) for x in sorted_indices])
        coords=self.frac_coords
        species = list(np.array(self.species)[sorted_indices])
        for i,site in enumerate(self):
            self[i]=species[i],coords[sorted_indices][i]

    def sort_sites(self):
        # First: Sort sites using the distance to the origin
        sorted_indices = np.array([np.linalg.norm(self.cart_coords[i]) for i in range(self.num_sites)]).argsort()
        self.sort_sites_using_list(sorted_indices)

        # Second: Sort again using the atomic number
        if len(self.types_of_specie) > 1:
            sorted_indices = np.array(self.atomic_numbers).argsort()
            self.sort_sites_using_list(sorted_indices)

    def sort_axes(self):
        """
        Sort the lattice vectors in decremental order of their size.
        'a' will be the shorest lattice vector
        'c' will be the longest one
        """

        sorted_indices = np.array(self.lattice.abc).argsort()
        #sorted_indices = np.array(self.lattice.abc).argsort()[::-1]
        self.modify_lattice(MLattice(self.lattice.matrix[sorted_indices]))
        frac_coords = self.frac_coords[:, sorted_indices]
        for i,site in enumerate(self):
            self[i]=site.specie,frac_coords[i]
  
    def align_with_axis(self, axis=0, round_decimals=14):
        lattice = MLattice(self.lattice.matrix)
        new_lattice=lattice.align_with_axis(axis=axis, round_decimals=round_decimals)
        #print(new_lattice)
        self.modify_lattice(MLattice(new_lattice))

    def align_with_plane(self, axis=2, round_decimals=14):
        lattice = MLattice(self.lattice.matrix)
        new_lattice=lattice.align_with_plane(axis=axis, round_decimals=round_decimals)
        self.modify_lattice(MLattice(new_lattice))

    def atoms_in_box(self):
        while min(self.frac_coords.flatten()) < 0.0 or max(self.frac_coords.flatten()) > 1.0:
            for i,site in enumerate(self):
                self.frac_coords[i] = site.specie,(frac_coords[i] + 1.0) % 1.0

    def canonical_form(self):
     #   self.sort_sites()
        self.sort_axes()
     #   print('----------')       #ok
     #   print(self.lattice.matrix)

        self.align_with_axis()
     #   print('----------')
     #   print(self.lattice.matrix)

        self.align_with_plane()
     #   print('----------')
     #   print(self.lattice.matrix)

        self.atoms_in_box()
        self.sort_sites()
                                             
class MMolecule(Molecule):

    @property
    def natom(self):
        assert self.composition.num_atoms - int(self.composition.num_atoms) < eps
        return int(self.composition.num_atoms)

    def canonical_form(self):
        self.relocate_to_cm()
        self.align_inertia_momenta()
        self.sort_sites()

    def align_inertia_momenta(self):
        I = self.inertia_matrix()
        eigval, eigvec = np.linalg.eig(I)
        eigvec = eigvec.T[eigval.argsort()[::-1]].T
        inveigvec = np.linalg.inv(eigvec)
        
       #  self.positions = np.dot(inveigvec, self.cart_coords.T).T
        coords=np.dot(inveigvec, self.cart_coords.T).T
        for i,site in enumerate(self):
            self[i]=site.specie,coords[i]

    def relocate_to_cm(self, list_of_sites=None):
        """
        Relocates the system of atoms to the center of mass
        a partial list of atoms can be used to compute
        the center, but all the atoms are moved to the
        computed center

        :param list_of_atoms: (list) List of atoms that will be considered for computing the center of mass
                              by default all atoms are included
        """
        cm = self.center_mass(list_of_sites)
        for i,site in enumerate(self):
            self[i]=site.specie,site.coords-cm

    def center_mass(self, list_of_sites=None):
        """
        Computes the center of mass (CM) of the XYZ object or
        a partial list of atoms. The default is to compute the
        CM of all the atoms in the object, if a list
        is enter only those in the list will be included for the CM
        Return the CM as a numpy array
        """
        if list_of_sites is None:
            list_of_sites = self.sites

        total_weight = 0.0
        center_of_mass = np.zeros(3)
        if self.natom == 0:
            return center_of_mass
        
        for site in self:
            if site in list_of_sites:
               wt = site.species_and_occu.weight
               center_of_mass += site.coords * wt
               total_weight += wt
        return center_of_mass / total_weight

    def inertia_matrix(self):
        im_xx = self.moment_of_inertia(0)
        im_yy = self.moment_of_inertia(1)
        im_zz = self.moment_of_inertia(2)
        im_xy = self.product_of_inertia(2)
        im_xz = self.product_of_inertia(1)
        im_yz = self.product_of_inertia(0)

        im = np.array([[im_xx, -im_xy, -im_xz], [-im_xy, im_yy, -im_yz], [-im_xz, -im_yz, im_zz]])
        return im

    def moment_of_inertia(self, axis):

        assert self.is_perfect
        I = 0
        for isite in self:
            I += isite.specie.data['Atomic mass'] * (sum(isite.coords ** 2) - isite.coords[axis] ** 2)
        return I

    def product_of_inertia(self, axis):

        assert self.is_perfect
        I = 0
        for isite in self:
            I += isite.specie.data['Atomic mass'] * (np.prod(isite.coords) / isite.coords[axis])
        return I

    @property
    def is_perfect(self):
       return  int(self.composition.num_atoms)==self.num_sites and min(self.composition.as_dict().values()) >= 1.0

    def sort_sites(self):

        # First: Sort sites using the distance to the origin
        sorted_indices = np.array([np.linalg.norm(self.cart_coords[i]) for i in range(self.num_sites)]).argsort()
        # print sorted_indices
        self.sort_sites_using_list(sorted_indices)

        # Second: Sort again using the atomic number
        if len(self.types_of_specie) > 1:
            sorted_indices = np.array(self.atomic_numbers).argsort()
            self.sort_sites_using_list(sorted_indices)

    def sort_sites_using_list(self, sorted_indices):
        sorted_indices = np.array([int(x) for x in sorted_indices])
        coords=self.cart_coords
        species = list(np.array(self.species)[sorted_indices])
        for i,site in enumerate(self):
            self[i]=species[i],coords[sorted_indices[i]]


if __name__=='__main__':

     
     #print('*'*30+'Molecule'+'*'*30)
     #mol= MMolecule.from_file('h2o.xyz')
     #print(mol)
     #print('-------canonical------------')
     #mol.canonical_form()
     #print(mol)
     print('*'*30+'Crystal'+'*'*30)
     st=MStructure.from_file('POSCAR')
     print(st)
     print('-------canonical------------')
     st.canonical_form()
     print(st)
     print('cell matrix')
     print(st.lattice.matrix)
     st.to(filename='p.vasp',fmt='poscar')
