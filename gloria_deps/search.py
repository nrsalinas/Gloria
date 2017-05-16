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

import data
import modelSel
import sys
import rescon
import hclust
from math import factorial, log, exp, pi
from itertools import combinations
from random import sample


def fieldPreSampler(observations, cohesion, debbug = False):
	preFields = {}
	bitRecord = []

	for ite in observations:
		medoid = ite
		cluster = []
		for o in observations:
			if data.getDist(o,medoid) < cohesion:
				cluster.append(o)
		if len(cluster) > 1:
			fieldie = data.HMRF(template=medoid)
			fieldie.emea(cluster)
			bits = fieldie.toBits()
			if bits not in bitRecord:
				#print bits
				bitRecord.append(bits)
				preFields[fieldie] = fieldie.pseudoLike(cluster)
	if debbug:
		assert len(preFields) > 0, "search.fieldPreSampler could not find a medoid."
		assert len(filter(lambda x: x.isNull() == True, preFields.keys())) < 1, "search.fieldPreSampler returned a null field as a medoid."
	return preFields

def reduceAreas(fieldDic, maxCombinations, numClusters):
	"""
	Reduces the list of probable fields given a maximum number of combinations
	allowed by the user.

	Arguments:

	- fieldDic: a dictionary with fields as keys and pseudolikelihoods as values.
	Typically, a dictionary returned by `fieldPreSampler` function.

	- maxCombinations: maximum number of combinations allowed.

	- numClusters: number of fields observations will be clustered to. In
	biogeographic context, the number of areas of endemism to be optimized.

	"""
	listLen = int()
	selAreas = []
	for k in xrange(len(fieldDic)):
		listLen = len(fieldDic) - k
		co = factorial(listLen) / (factorial(numClusters) * factorial(listLen - numClusters))
		if co <= maxCombinations:
			break

	selAreas = sorted(fieldDic, key = lambda x: fieldDic[x], reverse = True)
	selAreas = selAreas[:(listLen + 1)]
	return selAreas

def fieldOptim(observations, clusCohesion = 0.3, maxCycleIters = None, gammaParameter = 20, progress2stdout = False, pslikeFunc = 'raw', mixture = False, debbug = False):
	"""
	Main function wrapper. Execute model selection and optimize the state path
	and parameters of each HMRF.

	Arguments:

	- observations: list of data.Tile objects, usually modelled after actual
	taxa distributions.

	- clusCohesion (float): sets the dispersion or size of clusters during
	presampling.

	- maxCycleIters (int): maximum number of iterations (= maximum number of field
	combinations) allowed per model selection cycle.

	- gammaParameter (int): Potts model gamma parameter. Sets spatial correlation
	among neighbor cells in the field.
	"""

	selectedHypo = {}
	plics = []
	aics = []
	psdlks = []
	maxAreaNumber = 0
	result = rescon.Ensemble()

	# Test basic assumptions of function input
	if debbug:
		assert maxCycleIters is None or isinstance(maxCycleIters, int), "{0} is not a valid type for argument `maxCycleIters` in search.fieldOptim.".format(type(maxCycleIters))
		assert isinstance(gammaParameter , int) or isinstance(gammaParameter , float), "{0} is not a valid type for argument `gammaParameter` in search.fieldOptim.".format(type(gammaParameter))
		assert isinstance(clusCohesion , int) or isinstance(clusCohesion , float), "{0} is not a valid type for argument `clusCohesion` in search.fieldOptim.".format(type(clusCohesion))
		assert clusCohesion > 0 and clusCohesion < 1.0, "{0} is not a valid value for `clusCohesion` in search.fieldOptim function".format(clusCohesion)
		assert isinstance(observations , list), "`observations` argument in search.fieldOptim function is not a list."
		assert len(observations) > 0, "search.fieldOptim function called on no observations."
		assert len(observations) > 1, "search.fieldOptim function called on a single data.Tile."

	preFields = fieldPreSampler(observations, clusCohesion, debbug)
	if debbug:
		print "{0} seed fields".format(len(preFields))
		#for p in preFields:
		#	print p

	if progress2stdout:
		sys.stdout.write("\n")
		sys.stdout.flush()

	# Define maximum number of loops
	if len(observations) <= len(preFields):
		maxAreaNumber = len(observations) / 2
	else:
		maxAreaNumber = len(preFields)

	for k in xrange(1, (maxAreaNumber + 1)):
		if debbug:
			print "Assuming {0} components.".format(k)
		bestCluster = {}
		bestPseudolikelihood = -1e15
		plic = 0.0
		progStep = float()
		combCounter = 0

		if progress2stdout:
			backies = 36 + len(str(k)) + len(str(maxAreaNumber))
			sys.stdout.write("\b" * backies)
			sys.stdout.flush()
			sys.stdout.write("Cycle {0} of {1} posible =>  [          ]".format(k, maxAreaNumber))
			sys.stdout.flush()
			sys.stdout.write("\b" * 11)
			sys.stdout.flush()

		if maxCycleIters is None:
			newAreaSet = preFields.keys()
		else:
			newAreaSet = reduceAreas(preFields, maxCycleIters, k)

		progStep = (factorial(len(newAreaSet))/(factorial(k) * factorial(len(newAreaSet) - k))) / 10.0
		threshold = 0.0
		leftover = 10

		for comb in combinations(newAreaSet, k):
			if debbug:
				print "\tCombination",comb
			combCounter += 1
			if progress2stdout:
				if combCounter > threshold:
					sys.stdout.write(".")
					sys.stdout.flush()
					threshold += progStep
					leftover -= 1
			fields = []
			pseudolikelihood = 0.0
			#pseudolikelihoodMixture = 0.0
			for co in comb:
				fields.append(data.HMRF(template = co, gamma = gammaParameter))
			field2tiles = {x:[] for x in fields}
			for o in observations:
				bestDist = 1.0
				papa = None
				for f in fields:
					thisDist = data.getDist(o,f)
					if bestDist > thisDist and thisDist < clusCohesion:
						bestDist = thisDist
						papa = f
				if papa is not None:
					field2tiles[papa].append(o)

			solitaryFields = filter(lambda x : len(field2tiles[x]) < 2, field2tiles)
			if len(solitaryFields) > 0:
				# This optimization iteration does not contain all k clusters required
				continue

			if mixture:
				pseudolikelihood += modelSel.mixture_prob(field2tiles,pslikeFunc)
			else:
				for f in field2tiles:
					f.emea(field2tiles[f])
					if debbug:
						assert f.isNull() == False, "Expectation-Maximization algorithm resulted in a null field."
					pseudolikelihood += f.pseudoLike(field2tiles[f],pslikeFunc)

			if debbug:
				print "\tPseudolikelihood: {0}".format(pseudolikelihood)

			if bestPseudolikelihood < pseudolikelihood:
				bestPseudolikelihood = pseudolikelihood
				bestCluster = field2tiles

		# Should this loop be broken if no clusters were found?
		if bestCluster == {}:
			continue

		psdlks.append(bestPseudolikelihood)
		if mixture:
			plics.append(modelSel.plic_mixture(bestCluster, pslikeFunc))
		else:
			plics.append(modelSel.plic(bestCluster, pslikeFunc))

		aics.append(modelSel.aic(bestCluster)) ### Add AIC using mixtures

		if debbug:
			#print "\tPseudolikelihood = {0}".format(psdlks[-1])
			print "\tPLIC = {0}".format(plics[-1])

		if (len(plics) > 1 and plics[-1] < plics[-2]): #or (k == maxAreaNumber):
			if progress2stdout:
				if leftover > 0:
					sys.stdout.write("." * leftover)
					sys.stdout.write("]")
					sys.stdout.flush()

				sys.stdout.write("\n")
				sys.stdout.flush()
			if not debbug:
				break

		if k == maxAreaNumber and progress2stdout:
			if leftover > 0:
				sys.stdout.write("." * leftover)
				sys.stdout.write("]")
				sys.stdout.flush()

			sys.stdout.write("\n")
			sys.stdout.flush()

		else:
			result.field2taxa = bestCluster
			result.aic = aics[-1]
			result.plic = plics[-1]
			result.pseudolikelihood = psdlks[-1]
			result.components = len(bestCluster)

	# Retrieve noise
	for ob in observations:
		isNoise = True
		for fi in result.field2taxa:
			if ob in result.field2taxa[fi]:
				isNoise = False
				break
		if isNoise:
			result.noise.append(ob)

	return result

def hierarOptim(observations, maxDistance = 0.5, gamma = 20, pseudoLikeFunc = 'pmf', mixture = 'no', debbug = False):
	"""
	mixture: 'no' | 'complete' | 'simple'
	"""
	matrix = hclust.dist_matrix(observations, data.getDist)
	tree = hclust.dist_tree(matrix, 'UPGMA')
	#currentNodes = tree.node_filter(tree.root, maxdist = maxDistance)
	currentNodes = [tree.root]
	#print "Nodes_0: ",map(lambda x: sorted(tree.get_tips(x)), currentNodes)
	plicValues = []
	pseValues = []
	aicValues = []
	areaCounter = []
	node_maps = []
	result = rescon.Ensemble()
	nextGen = True
	prevGenNodes = len(currentNodes)

	if currentNodes == []:
		if debbug:
			print "Clustering step failed!"
		nextGen = False

	while nextGen:
		if debbug:
			print "\nElements:",len(currentNodes)
		fieldTaxDic = {}
		gotNullField = False
		nullFieldCount = 0

		for nod in currentNodes:
			field = hclust.clade2field(tree, nod, observations, gamma_in = gamma)
			#if debbug:
			#	print id(field)
			taxa = [observations[x] for x in tree.get_tips(nod)]
			#if debbug:
				#print "Indexes:",sorted(tree.get_tips(nod))
				#print "Cluster size:",len(tree.get_tips(nod))
			fieldTaxDic[field] = taxa

			if field.isNull():
				#if debbug:
				#	print "Field is null"
				field.setMeans(0.5,0.5)
				field.setStaDevs(0.01,0.00999)
				gotNullField = True
				nullFieldCount += 1

		net_clusters = len(currentNodes) - nullFieldCount
		if debbug:
			print "net_clusters:",net_clusters

		if net_clusters > 0:# gotNullField == False:

			areaCounter.append(len(currentNodes))
			if mixture == 'complete':
				plicValues.append(modelSel.plic_mixture(fieldTaxDic, pslikeplic = pseudoLikeFunc))
			elif mixture == 'no':
				plicValues.append(modelSel.plic(fieldTaxDic, plikefunc = pseudoLikeFunc, ms_debbug = debbug))
			elif mixture == 'simple':
				plicValues.append(modelSel.plic_weigthed(fieldTaxDic, pslikeplic = pseudoLikeFunc))
			if debbug:
				print "PLIC: ",plicValues[-1]
			aicValues.append(modelSel.aic(fieldTaxDic))
			if mixture == 'complete':
				pseudolikelihood = modelSel.mixture_prob(fieldTaxDic, pslikeFunc = pseudoLikeFunc)
			elif mixture == 'no':
				pseudolikelihood = sum([fi.pseudoLike(fieldTaxDic[fi], distr_form = pseudoLikeFunc) for fi in fieldTaxDic])
			elif mixture == 'simple':
				pseudolikelihood = modelSel.weighted_prob(fieldTaxDic, pslikeFunc = pseudoLikeFunc)
			pseValues.append(pseudolikelihood)
			node_maps.append(currentNodes)

		if len(currentNodes) >= (len(observations)/ 2):
			nextGen = False

		currentNodes = tree.clade_budding(currentNodes, maxDistance)

		if len(currentNodes) == prevGenNodes:
			if debbug:
				print "Clade budding failed"
			nextGen = False
		prevGenNodes = len(currentNodes)

	# Reversely find the first plic value peak
	alphaInd = None
	prevPlic = -1e44
	for idx, pl in reversed(list(enumerate(plicValues))):
		if pl > prevPlic:
			prevPlic = pl
			alphaInd = idx
		else:
			break
	#alphaInd = plicValues.index(max(plicValues))
	if debbug:
		print "plicValues:", plicValues
		print "alphaInd:", alphaInd
		#print "node_maps:", node_maps

	fieldTaxDic = {}
	for nod in node_maps[alphaInd]:
		field = hclust.clade2field(tree, nod, observations, gamma_in = gamma)
		if not field.isNull():
			result.field2taxa[field] = [observations[x] for x in tree.get_tips(nod)]
	result.aic = aicValues[alphaInd]
	result.plic = plicValues[alphaInd]
	result.pseudolikelihood = pseValues[alphaInd]
	result.components = len(result.field2taxa)

	return result
