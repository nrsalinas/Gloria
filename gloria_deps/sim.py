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
#from data import Tile
import random
from data import Tile

def getFake(case, clus=2, num=3, datRatio = 1, inun=0.95, exun=0.95, offun=0, noise=0, height=10, width=10, geometry="square"):
	"""
	Returns simulated geographic distributions as a list of data.Tile objects.
	Distributions are based on simulations presented by Casagranda et al. 2012,
	Cladistics 28: 645--654.

	Arguments:

	- case (int, 0-3): kind of simulation. Could be 0 (non-conflictive areas), 1
	(nested areas), 2 (overlapping areas), and 3 (disjoint areas). See Casagranda
	et al. for a detailed explanation.

	- clus (int, 1-3): number of clusters to simulate.

	- num (int or list of ints): number of datasets per cluster to simulate. If
	an integer is parsed (2-any), all clusters will have the same number of
	datasets. If a list of integers is parsed, clusters will have different sizes.
	Thus, parsing `[2,6,12]` will simulated 2 distributions for cluster 0, 6 for
	cluster 1, and 12 for cluster 2. The list length should be equal to the number
	of clusters set with argument `clus`.

	- inun (float, 0.0-1.0): degree of internal uncertainty within distributions
	(0.0 == maximum uncertainty, 1.0 == maximum certainty).

	- exun (float, 0.0-1.0): degree of internal uncertainty within distributions
	(0.0 == maximum uncertainty, 1.0 == maximum certainty).

	- offun (int, 1-3): degree of overlapping uncertainty among distributions
	belonging to the same cluster.

	- noise (int): number of non-clustering random distributions.

	- height (int): number of rows in the distribution lattice.

	- width (int): number of columns in the distribution lattice.

	- geometry (string, "square" or "hexagon"): geometric shape of cells making
	up the lattice.

	"""
	outList = []
	#listShape = [range((num * x),(num * (x + 1))) for x in xrange(clus)] # how outList indexes should be clustered
	clusterCount = 1

	xScale = {0.3:0, 0.5:0, 0.7:0}
	yScale = {0.3:0, 0.5:0, 0.7:0}
	for factor in xScale:
		xScale[factor] = int(factor * width) - 1
	for factor in yScale:
		yScale[factor] = int(factor * height) - 1
	x0 = 0
	xf = 0
	y0 = 0
	yf = 0
	noiseSize = 0
	neighborhood = int()
	if type(num) == int:
		obs = [num for x in xrange(clus)]
	elif type(num) == list:
		obs = num
	if case == 0: ### Non-problematic areas
		for c,o in zip(xrange(clus), obs):
			for n in xrange(o):
				matrix = [[0 for x in xrange(height)] for x in xrange(width)]
				xoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				yoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				neighborhood = int(yScale[0.5] * 0.25)
				if c == 0:
					x0 = 0 + xoffset
					xf = xScale[0.5] + xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] + yoffset
				elif c == 1:
					x0 = xScale[0.5] + 1 - xoffset
					xf = (width - 1) - xoffset
					y0 = yScale[0.5] + 1 - yoffset
					yf = (height - 1) - yoffset
				elif c == 2:
					x0 = 0 + xoffset
					xf = xScale[0.5] + xoffset
					y0 = yScale[0.5] + 1 - yoffset
					yf = (height - 1) - yoffset
				elif c == 3:
					x0 = xScale[0.5] + 1 - xoffset
					xf = (width - 1) - xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] + yoffset
				for ir,r in enumerate(matrix):
					for ic,cell in enumerate(r):
						myrand = numpy.random.random_sample()
						if (y0 <= ir <= yf) and (x0 <= ic <= xf):
							if myrand < inun:
								matrix[ir][ic] = 1
						elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
							if myrand > exun:
								matrix[ir][ic] = 1
				tilita = Tile(ingrid = matrix, cellType = geometry)
				outList.append(tilita)



	elif case == 1: ### Nested areas
		for c,o in zip(xrange(clus), obs):
			for n in xrange(o):
				matrix = [[0 for x in xrange(height)] for x in xrange(width)]
				xoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				yoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				if c == 0:
					x0 = 0 + xoffset
					xf = (width - 1) - xScale[0.3] + xoffset
					y0 = 0 + yoffset
					yf = (height - 1) - yScale[0.3] + yoffset
					neighborhood = int((height - yScale[0.3]) * 0.25)
				elif c == 1:
					x0 = 0 + xoffset
					xf = xScale[0.5] + xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] + yoffset
					neighborhood = int(xScale[0.5] * 0.25)
				elif c == 2:
					x0 = 0 + xoffset
					xf = xScale[0.3] + xoffset
					y0 = 0 + yoffset
					yf = yScale[0.3] + yoffset
					neighborhood = int(xScale[0.3] * 0.25)
				for ir,r in enumerate(matrix):
					for ic,cell in enumerate(r):
						myrand = numpy.random.random_sample()
						if (y0 <= ir <= yf) and (x0 <= ic <= xf):
							if myrand < inun:
								matrix[ir][ic] = 1
						elif (ir <= (yf + neighborhood)) and (ic <= (xf + neighborhood)):
							if myrand > exun:
								matrix[ir][ic] = 1
				tilita = Tile(ingrid = matrix, cellType = geometry)
				outList.append(tilita)
	elif case == 2: ### Overlapping areas
		for c,o in zip(xrange(clus), obs):
			for n in xrange(o):
				matrix = [[0 for x in xrange(height)] for x in xrange(width)]
				xoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				yoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				neighborhood = int(yScale[0.5] * 0.25)
				if c == 0:
					x0 = 0 + xoffset
					xf = xScale[0.5] + xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] + yoffset
				elif c == 1:
					x0 = xScale[0.3] + xoffset
					xf = xScale[0.5] + xScale[0.3] + xoffset
					y0 = xScale[0.3] + yoffset
					yf = yScale[0.5] + xScale[0.3] + yoffset
				elif c == 2:
					x0 = xScale[0.5] + 1 - xoffset
					xf = (width - 1) - xoffset
					y0 = yScale[0.5] + 1 - yoffset
					yf = (height - 1) - yoffset
				for ir,r in enumerate(matrix):
					for ic,cell in enumerate(r):
						myrand = numpy.random.random_sample()
						if (y0 <= ir <= yf) and (x0 <= ic <= xf):
							if myrand < inun:
								matrix[ir][ic] = 1
						elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
							if myrand > exun:
								matrix[ir][ic] = 1
				tilita = Tile(ingrid = matrix, cellType = geometry)
				outList.append(tilita)
	elif case == 3: ### Disjoint areas
		for c,o in zip(xrange(clus), obs):
			for n in xrange(o):
				matrix = [[0 for x in xrange(height)] for x in xrange(width)]
				xoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				yoffset = int(offun * numpy.random.normal(1,0.2,1)[0])
				if c == 0:
					neighborhood = int(yScale[0.5] * 0.25)
					x0 = 0 + xoffset
					xf = xScale[0.5] - 1 + xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] - 1 + yoffset
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
					x0 = xScale[0.5] + 2 - xoffset
					xf = (width - 1) - xoffset
					y0 = yScale[0.5] + 2 - yoffset
					yf = (height - 1) - yoffset
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
				elif c == 1:
					neighborhood = int(yScale[0.5] * 0.25)
					x0 = 0 + xoffset
					xf = xScale[0.5] - 1 + xoffset
					y0 = yScale[0.5] + 2 - yoffset
					yf = (height - 1) - yoffset
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
					x0 = xScale[0.5] + 2 - xoffset
					xf = (width - 1) - xoffset
					y0 = 0 + yoffset
					yf = yScale[0.5] - 1 + yoffset
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
				elif c == 2:
					neighborhood = int(yScale[0.5] * 0.25)
					x0 = xScale[0.3] + xoffset
					xf = x0 + xScale[0.5]
					y0 = 0 + yoffset
					yf = y0 + yScale[0.3]
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
					x0 = xScale[0.3] - xoffset
					xf = x0 + xScale[0.5]
					y0 = yScale[0.5] + yScale[0.3] - yoffset
					yf = (height - 1) - yoffset
					for ir,r in enumerate(matrix):
						for ic,cell in enumerate(r):
							myrand = numpy.random.random_sample()
							if (y0 <= ir <= yf) and (x0 <= ic <= xf):
								if myrand < inun:
									matrix[ir][ic] = 1
							elif ((y0 - neighborhood) <= ir <= (yf + neighborhood)) and ((x0 - neighborhood) <= ic <= (xf + neighborhood)):
								if myrand > exun:
									matrix[ir][ic] = 1
				tilita = Tile(ingrid = matrix, cellType = geometry)
				outList.append(tilita)

	ind = range(height) + range(height, (height + width) )
	sel = random.sample(ind,noise)
	for n in sel:
		if n < height:
			matrix = [[0 for x in xrange(height)] for x in xrange(width)]
			matrix[n] = [1 for x in xrange(width)]
			tilita = Tile(ingrid = matrix, cellType = geometry)
			outList.append(tilita)
		else:
			c = n - height
			matrix = [[0 for x in xrange(height)] for x in xrange(width)]
			for ir,row in enumerate(matrix):
				for ic,col in enumerate(row):
					if ic == c:
						matrix[ir][ic] = 1
			tilita = Tile(ingrid = matrix, cellType = geometry)
			outList.append(tilita)

	return outList

def getBigDataset(clus = 30, num = 2, inun=0.95, exun=0.95, noise = 0):
	"""
	Returns simulated geographic distributions as a list of data.Tile objects.
	Lattices are 30x30 cell grids, effective distributions are 5x5.

	Arguments:

	clus (int, 1-36): Number of clusters to simulate.

	num (int, 2-any): Number of observations to simulate per cluster.

	inun (float, 0.0-1.0): degree of internal uncertainty within distributions
	(0.0 == maximum uncertainty, 1.0 == maximum certainty).

	exun (float, 0.0-1.0): degree of internal uncertainty within distributions
	(0.0 == maximum uncertainty, 1.0 == maximum certainty).

	"""
	universeSize = 30
	areasSize = 5
	cluscount = 0
	breakall = False
	thisyoff, thisxoff = 0, 0
	outList = []
	internalCellsOff = int((1-inun) * areasSize**2)
	externalCellsOn = int((1-exun) * areasSize**2)
	for yoff in xrange(6):
		if breakall:
			break
		for xoff in xrange(6):
			cluscount += 1
			if cluscount > clus:
				breakall = True
				break
			thisyoff, thisxoff = yoff, xoff
			for x in xrange(num):
				listmp = [[0 for x in xrange(universeSize)] for x in xrange(universeSize)]
				for ir in xrange((thisyoff * areasSize), ((thisyoff + 1) * areasSize)):
					for ic in xrange((thisxoff * areasSize), ((thisxoff + 1) * areasSize)):
						listmp[ir][ic] = 1
				yzero = random.sample(range((thisyoff * areasSize), ((thisyoff + 1) * areasSize)),internalCellsOff)
				xzero = random.sample(range((thisxoff * areasSize), ((thisxoff + 1) * areasSize)),internalCellsOff)

				yones, xones = [], []

				while len(yones) < externalCellsOn:
					cand = random.sample(range(((thisyoff * areasSize) - 1), (((thisyoff + 1) * areasSize) + 1)), 1)[0]
					if cand >= 0 and cand < 30:
						yones.append(cand)

				while len(xones) < externalCellsOn:
					cand = random.sample(range(((thisxoff * areasSize) - 1), (((thisxoff + 1) * areasSize) + 1)), 1)[0]
					if cand >= 0 and cand < 30:
						ti = len(xones)
						if (thisyoff * areasSize) <= yones[ti] and (thisxoff + 1) * areasSize > yones[ti]:
							if cand >= ((thisxoff + 1) * areasSize) or cand < (thisxoff * areasSize):
								xones.append(cand)
						else:
							xones.append(cand)
				for y in xrange(internalCellsOff):
					listmp[yzero[y]][xzero[y]] = 0
				for y in xrange(externalCellsOn):
					listmp[yones[y]][xones[y]] = 1
				outList.append(Tile(ingrid=listmp,cellType='square'))

	if noise > 0:
		indexes = range(2 * universeSize)
		rands = random.sample(indexes,noise)
		for mynoi in rands:
			listmp = [[0 for x in xrange(universeSize)] for x in xrange(universeSize)]
			if mynoi >= universeSize:
				for i in range(universeSize):
					listmp[i][(mynoi - universeSize)] = 1
			else:
				for i in range(universeSize):
					listmp[mynoi][i] = 1
			outList.append(Tile(ingrid=listmp,cellType='square'))
	return outList



#datasets = getMat(5,2)
