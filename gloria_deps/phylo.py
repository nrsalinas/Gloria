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

class Tree(object):
	def __init__(self, num_leaves):
		# Init lists in tree with max number of nodes given leaes
		self.total_nodes = 2 * num_leaves - 1
		# Reference lists
		self.nodes = [x for x in xrange(self.total_nodes)]
		self.parents = [-1 for x in xrange(self.total_nodes)]
		self.old_children = [-1 for x in xrange(self.total_nodes)]
		self.siblings = [-1 for x in xrange(self.total_nodes)]
		self.branches = [0.0 for x in xrange(self.total_nodes)]
		self.distance_to_tips = [0.0 for x in xrange(self.total_nodes)]
		self.offspring = [0 for x in xrange(self.total_nodes)]
		# Key variables. Refer to indexes in self.nodes
		self.root = -2
		self.lastLeaf = num_leaves - 1
		self.lastNode = self.lastLeaf
		# For scipy plotting
		self.linking_matrix = []

	def children_gen(self, node):
		node = self.old_children[node]
		if node != -1:
			yield node
			node = self.siblings[node]
			while node != -1:
				yield node
				node = self.siblings[node]

	def gen(self, node):
		if self.siblings[node] != -1:
			for sib in self.gen(self.siblings[node]):
				yield sib
		if self.old_children[node] != -1:
			for lc in self.gen(self.old_children[node]):
				yield lc
		yield node

	def check_integrity(self):
		out = True
		# Check all nodes (but root) have parents
		if self.parents[self.root] != -1:
			out = False
		dummy = [x for x in self.nodes]
		dummy.remove(self.root)
		orphans = filter(lambda x : self.parents[x] == -1, dummy)
		if orphans != []:
			out = False

		# Check all nodes (but leaves) have old_children
		leaves = self.nodes[:(self.lastLeaf + 1)]
		inter_nodes = self.nodes[(self.lastLeaf + 1):]
		early_parenting = filter(lambda x : self.old_children[x] != -1, leaves)
		if early_parenting != []:
			out = False
		bachelors = filter(lambda x : self.old_children[x] == -1, inter_nodes)
		if bachelors != []:
		 	out = False

		# Check all nodes can be iterated from root.
		iternodes = [n for n in self.gen(self.root)]
		if sorted(self.nodes) != sorted(iternodes):
			out = False

		# Check all nodes are older than descendants
		for n in self.gen(self.root):
			if self.offspring[n] > 1:
				for ch in self.children_gen(n):
					if self.distance_to_tips[n] - self.distance_to_tips[ch] <= 0:
						out = False
		return out

	def get_tips(self, node):
		tips = []
		for child in self.children_gen(node):
			if self.old_children[child] == -1:
				tips.append(child)
			else:
				tips += self.get_tips(child)
		if tips == []:
			tips = [node]
		return tips

	def node_filter(self, node, maxdist = 1.0):
		"""
		Returns a list of nodes which is distance to the tip is not greater than
		`maxdist`. The list do not include descendants from nodes already appended.
		"""
		filteredNodes = []
		if self.distance_to_tips[node] > maxdist:
			if self.offspring[node] > 0:
				for child in self.children_gen(node):
					filteredNodes += self.node_filter(child, maxdist)
		else:
			if self.offspring[node] > 0:
				filteredNodes.append(node)

		return filteredNodes


	def clade_budding(self, nodes, maxdist = 0.5):
		out = []
		to_rm = []
		newnodes = []
		innodes = len(nodes)

		nodes.sort(reverse = True, key = lambda x: self.distance_to_tips[x])
		for ind, n in enumerate(nodes):
			descendants = [self.offspring[x] for x in self.children_gen(n)]
			non_leaves = filter(lambda x: x > 0, descendants)

			if self.distance_to_tips[n] < maxdist and \
						len(non_leaves) == len(descendants):
				to_rm.append(ind)
				break

			elif self.distance_to_tips[n] > maxdist and \
						len(non_leaves) > 0:
				# Probably this will cause to return a node count greater
				# than (one + the current size)
				to_rm.append(ind)
				if len(non_leaves) > 1:
					break
		out = [n for n in nodes]
		for r in to_rm:
			bye = out.pop(r)
			out += filter(lambda c: self.offspring[c] > 0, \
				self.children_gen(bye))

		return out
