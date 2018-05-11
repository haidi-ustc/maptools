#!/usr/bin/env python
from   __future__ import division, unicode_literals, print_function
import math
import numpy as np
import scipy.spatial
import itertools
from   collections import OrderedDict
from   pymatgen.core.periodic_table import Element
from   pymatgen import Molecule
from   maptool.maptool_configuration import MMolecule,MStructure
from   maptool.mathematics import integral_gaussian
from   maptool import mpt_log

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2018"
__version__ = "0.17.3"
__email__ = "gufranco@mail.wvu.edu"
__status__ = "Development"
__date__ = "March 7, 2017"
__maintainer__="haidi Wang"

class MoleculeAnalysis:
    def __init__(self, molecule):
        """
        Cluster Analysis provides routines to compute Molecule Analysis
        for finite systems such as molecules and clusters.

        :param molecule: (pymatgen.Molecule) A Pymatgen Molecule object

        """
        cm=molecule.center_of_mass
        coords=molecule.cart_coords-cm
        self.molecule=MMolecule(molecule.species,coords)
        # This is the maximal distance between any two atoms
        self.max_distance = np.max(self.molecule.distance_matrix)
        self._distances = None
        self._all_distances = None
        self._pairs = None


    def close_distances(self):
        """
        Computes the closest distances for all the atoms

        :return: (tuple) Return a bond's dictionary and distance's list
        """
        if self._pairs is None or self._distances is None:
           dm = self.molecule.distance_matrix
           dm += np.eye(len(dm)) * max(dm.flatten())
           pairs_dict = {}
           distances_list = []
           for i in range(self.molecule.natom):
               index = dm[:, i].argmin()
               pairs_dict[str(i)] = [index]
               distances_list.append(dm[index, i])
           self._pairs = pairs_dict
           self._distances = distances_list

        return self._pairs, self._distances

    def all_distances(self):

        if self._all_distances is None:
            ret = {}
            dist_matrix = self.molecule.distance_matrix
            for i, j in itertools.combinations_with_replacement(range(self.molecule.natom), 2):
                pair = (i, j)
                ret[pair] = dist_matrix[i, j]
            self._all_distances = ret

        return self._all_distances


    def distance_matrix(self):
        """
        Compute the distance matrix, the distance between any two particles
        in the entire molecule.
        For distance matrices related separated by species, use instead
        distance_matrix_by_species.

        :return:
        """
        return self.molecule.distance_matrix

    def distance_matrix_by_species(self, specie1, specie2):

#        assert self.molecule.is_perfect
        group1 = np.array(self.molecule.species) == Element(specie1)
        group2 = np.array(self.molecule.species) == Element(specie2)
        return scipy.spatial.distance_matrix(self.molecule.cart_coords[group1], self.molecule.cart_coords[group2])

    def all_distances_by_species(self):

        ret = OrderedDict()
        for i, j in itertools.combinations_with_replacement(self.molecule.types_of_specie, 2):
            pair = tuple(np.sort([i, j]))
            dm = self.distance_matrix_by_species(i, j)
            if dm.shape[0]  ==1 and dm.shape[1] !=1:
               distances=dm[0]
            elif dm.shape[0]!=1 and dm.shape[1] ==1:
               distances=dm.T[0]
            elif dm.shape[0]==1 and dm.shape[1] ==1:
               distances = dm[0]
            else:
               distances = np.triu(dm, 1).flatten()
               distances = np.sort(distances[distances > 0.0])
            ret[pair] = distances
        return ret

    def discrete_radial_distribution_function(self, delta=0.01, sigma=0.01, integrated=True):
        dist_spec_bk = self.all_distances_by_species()
        dist_spec=dist_spec_bk.copy()
        for key in dist_spec_bk.keys():
            if key[0]==key[1] and dist_spec_bk[key][0]==0.0:
               del dist_spec[key]
        discrete_rdf = OrderedDict()
        nbins = int((self.max_distance + 8 * sigma) / delta)
        if nbins > 10000:
            nbins = 10000
        discrete_rdf_x = np.arange(0, nbins * delta, delta)
        for spec_pair in dist_spec:
            discrete_rdf[spec_pair] = np.zeros(nbins)
            for Rij in dist_spec[spec_pair]:
                # Integrating for bin from - 8*sigma to +8*sigma centered on Rij
                # Values outside this range are negligible
                imin = int(max(0, (Rij - 8 * sigma) / delta))
                imax = int(min(len(discrete_rdf_x), (Rij + 8 * sigma) / delta))
                mpt_log.debug('Computing for distance %7.3f for indices between %d and %d' % (Rij, imin, imax))
                for i in range(imin, imax):
                    x = discrete_rdf_x[i]
                    if not integrated:
                        discrete_rdf[spec_pair][i] += np.exp(-((x - Rij) ** 2) / (2 * sigma * sigma)) / (
                            4 * math.pi * Rij * Rij)
                    else:
                        discrete_rdf[spec_pair][i] += integral_gaussian(x, x + delta, Rij, sigma) / (
                            4 * math.pi * Rij * Rij)

        return discrete_rdf_x, discrete_rdf

    def fp_oganov(self, delta=0.01, sigma=0.01):
        struc_dist_x, struc_dist = self.discrete_radial_distribution_function(delta=delta, sigma=sigma)
        fp_oganov=OrderedDict()
        key_order=sorted([tuple(sorted(list(x))) for x in struc_dist.keys()])
        for spec_pair in key_order:
            new_key=(spec_pair[0].number,spec_pair[1].number)
            fp_oganov[new_key]=struc_dist[spec_pair]
        return struc_dist_x, fp_oganov


class ClusterMatch:
    def __init__(self, molecule1, molecule2):
        """
        Given two molecules, ClusterMatch can compute a one-to-one
        association between atoms of both molecules, trying to minimize
        the displacements needed to move one molecule into another.

        :param molecule1:
        :param molecule2:
        """

        assert isinstance(molecule1,MMolecule)
        assert isinstance(molecule2,MMolecule)

        self.molecule1 = molecule1.copy()
        self.molecule2 = molecule2.copy()

        self.molecule1.canonical_form()
    
        self.molecule2.canonical_form()

        assert self.molecule1.species == self.molecule2.species
        self._matched=False

    def match(self):

        species = self.molecule1.types_of_specie
        species = list(np.array(species)[np.array(self.molecule1.atomic_numbers).argsort()])

        permutation = np.zeros(len(self.molecule1.cart_coords), dtype=int)
        num = 0
        for ispecie in species:

            pos1 = self.molecule1.cart_coords[np.array(self.molecule1.types_of_specie) == ispecie]
            pos2 = self.molecule2.cart_coords[np.array(self.molecule2.types_of_specie) == ispecie]
            dm = scipy.spatial.distance_matrix(pos1, pos2)
            match_list = np.zeros(len(pos1))
            maxdis = np.max(dm.flatten())
            pcm_log.debug(dm)
            for i in range(len(pos1)):
                match = dm[i, :].argmin()   # return the index of minimum value
                # print " %3s %3d %9.3f %9.3f %9.3f" % (ispecie, match, np.linalg.norm(pos1[i]),
                # np.linalg.norm(pos1[match]), dm[i,match])
                match_list[i] = match
                dm[:, match] = maxdis + 1

            match_list += num
            permutation[num:num + len(pos1)] = match_list
            num += len(pos1)

        # print permutation
        self.molecule2.sort_sites_using_list(permutation)
        self._matched=True
      
    def match_result(self):
        if self._matched:
           return self.molecule1,self.molecule2
        else:
           self.match()
           return self.molecule1,self.molecule2 
