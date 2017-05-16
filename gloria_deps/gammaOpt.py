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

from random import sample, choice, randrange
from data import getDist, HMRF
from math import exp, log

def spaCorr(obs, iterations = 100):
	distr = {}
	for x in xrange(iterations):
		thisObs = choice(obs)
		thisCol = randrange(0, obs[0].cols)
		thisRow = randrange(0, obs[0].rows)
		thisSymbol = thisObs[thisRow, thisCol]
		neighbors = thisObs.getNeighs(thisRow, thisCol)
		same = 0.0
		for neir,neic in neighbors:
			if thisSymbol == thisObs[neir, neic]:
				same += 1.0
		same /= float(len(neighbors))
		same = round(same, 2)
		if same in distr:
			distr[same] += 1
		else:
			distr[same] = 1
	distr = {key : (distr[key] / float(iterations)) for key in distr}
	return distr

def ufunc(gamma, homogeneity):
	if homogeneity >= 0 and homogeneity <= 1 and gamma > 0:
		return exp(homogeneity * gamma) / (exp(homogeneity * gamma) + exp((1.0 - homogeneity) * gamma))

def lnFit(homoDic, gammaList = None):
	mygamma = 20
	myerror = 1e10

	if gammaList is None:
		gammaList = range(1,31)
	distances = {ga: 0.0 for ga in gammaList}

	for gamma in gammaList:
		for homogeneity in homoDic:
			expect = ufunc(gamma, homogeneity)
			distances[gamma] += (log(expect) - log(homoDic[homogeneity])) ** 2
		#print distances[gamma]

	for g in reversed(sorted(distances)):
		if myerror > distances[g]:
			myerror = distances[g]
			mygamma = g

	return mygamma
