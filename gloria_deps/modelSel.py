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

import numpy
from math import factorial, log, exp, pi

def aic(dic):
	"""
	Receives a dictionary of Tile objects, typically the dic attribute of a
	Cluster object. Keys are a Markov Random Fields, values are list of cluster
	elements.
	"""
	aic = 0.0
	numpars = float(len(dic)) * 5.0 # mean and standard deviation per each state plus a parameter for the prior energy function
	for key in dic:
		assert len(dic[key]) > 0, "Cluster has no elements (modelSel.aic)."
		aic += key.pseudoLike(dic[key])
	aic = ((-2.0) * aic) + (2.0 * numpars)
	return aic

def plic(dic , plikefunc, mixture = False, ms_debbug = False):
	"""
	Receives a dictionary of Tile objects, typically the dic attribute of a
	Cluster object. Keys are Markov Random Fields, values are list of cluster
	elements.
	"""
	nkeys = len(dic)
	mylambda = 3.0
	plicv = 0.0
	numpars = float(len(dic)) * 5.0 # mean and standard deviation per each state plus a parameter for the prior energy function
	N = 0.0
	allObs = reduce(lambda x,y: x+y, dic.values())
	for key in dic:
		if ms_debbug:
			print "Field {0}, observations: {1},".format(id(key), len(dic[key])),
		assert len(dic[key]) > 0, "Cluster has no elements (modelSel.plic)."
		if ms_debbug and key.isNull():
			print "null field, ",
		if mixture:
			thisPL = key.pseudoLike(allObs, plikefunc)
		else:
			thisPL = key.pseudoLike(dic[key], plikefunc)
		if ms_debbug:
			print "pseudolikelihood: {0}".format(thisPL)
		plicv += thisPL
		N += float(len(dic[key]) * key.rows * key.cols)

	if len(dic) == 0:
		plicv = 0.0
	else:
		plicv = (2 * plicv) - (log(N) * numpars)
	return plicv


def weighted_prob(dic, pslikeFunc = 'raw'):
	nelems = len(reduce(lambda x,y: x+y, dic.values()))
	pseudo = 0.0
	for field in dic:
		weight = len(dic[field]) / float(nelems)
		this_prob = log(weight) + field.pseudoLike(dic[field], distr_form = pslikeFunc)
		#print "this logprob:",this_prob
		pseudo += this_prob
	return pseudo


def mixture_prob(dic, pslikeFunc = 'raw'):
	nelems = len(reduce(lambda x,y: x+y, dic.values()))
	pseudo = 0.0
	for field in dic:
		weight = len(dic[field]) / float(nelems)
		for obs in dic[field]:
			this_prob = log(weight) + field.pseudoLike([obs], distr_form = pslikeFunc)
			pseudo += exp(this_prob)
		try:
			pseudo = log(pseudo)
		except:
			pseudo = -745.0
	return pseudo


def plic_mixture(dic, pslikeplic = 'raw'):
	logprob = mixture_prob(dic, pslikeplic)
	N = float(sum([(len(dic[key]) * key.rows * key.cols) for key in dic]))
	numpars = float(len(dic)) * 5.0 # mean and standard deviation per each state plus a parameter for the prior energy function
	plic = (2 * logprob) - (log(N) * numpars)
	return plic

def plic_weigthed(dic, pslikeplic = 'raw'):
	logprob = weighted_prob(dic, pslikeplic)
	N = float(sum([(len(dic[key]) * key.rows * key.cols) for key in dic]))
	numpars = float(len(dic)) * 5.0 # mean and standard deviation per each state plus a parameter for the prior energy function
	plic = (2 * logprob) - (log(N) * numpars)
	return plic
