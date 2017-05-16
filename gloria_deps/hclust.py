###############################################################################
#
#	Copyright 2016-2017 Nelson R. Salinas
#
#
#	This file is part of Gloria.
#
#   Gloria is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#
#	Gloria is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with Gloria.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

import fastcluster
import phylo
from data import getDist, HMRF

def dist_tree(distance_matrix, link_alg = 'UPGMA'):
	if link_alg == 'UPGMA':
		z = fastcluster.linkage(distance_matrix, method='average')

	if len(z) > 0:
		tree = phylo.Tree(len(distance_matrix))
		tree.linking_matrix = z
		for row in z:
			left = int(row[0])
			right = int(row[1])
			tree.lastNode += 1
			tree.parents[left] = tree.lastNode
			tree.parents[right] = tree.lastNode
			tree.old_children[tree.lastNode] = left
			tree.siblings[left] = right
			tree.offspring[tree.lastNode] = 2
			tree.distance_to_tips[tree.lastNode] = (tree.distance_to_tips[left] \
				+ tree.distance_to_tips[right] + row[2]) / 2
			tree.branches[left] = tree.distance_to_tips[tree.lastNode] - \
				tree.distance_to_tips[left]
			tree.branches[right] = tree.distance_to_tips[tree.lastNode] - \
				tree.distance_to_tips[right]
		tree.root = tree.lastNode
	return tree

def dist_matrix(observations, distance_function):
	"""
	Returns a distance matrix from a list of Tile objects. Argument 0 is a list
	of objects, argument 1 is a function to measure their pairwise distance.
	"""
	mat = [[x for x in xrange(len(observations))]
			for y in xrange(len(observations))]
	for row in xrange(len(observations)):
		for col in xrange(len(observations)):
			mat[row][col] = distance_function(observations[row],
									observations[col])
	return mat

def clade2field(tree, node, observations, gamma_in):
	leaves = tree.get_tips(node)
	obsSubset = [observations[l] for l in leaves]
	field = HMRF(ingrid=[[0 for x in xrange(observations[0].cols)] \
			for y in xrange(observations[0].rows)], gamma = gamma_in)
			#, mean0 = 0.4, mean1 = 0.5)
	field.emea(obsSubset)
	return field
