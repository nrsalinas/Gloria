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

from libc.math cimport exp, log, M_PI
import random
import numpy as np

cdef double getDistance(Tile tileA, Tile tileB) except -1:
	cdef:
		int ir, ic
		double dist, shared = 0.0, countA = 0.0, countB = 0.0
	for ir in xrange(tileA.rows):
		for ic in xrange(tileA.cols):
			if tileA.mvsymbols[ir,ic] == 1.0:
				countA += 1.0
			if tileB.mvsymbols[ir,ic] == 1.0:
				countB += 1.0
			if tileA.mvsymbols[ir,ic] == 1.0 and tileB.mvsymbols[ir,ic] == 1.0:
				shared += 1.0
	dist = (1.0 - (0.5 * ((float(shared) / float(countA)) + (float(shared) / float(countB)))))
	return dist

cdef double euclidean(Tile tileA, Tile tileB) except -1:
	cdef:
		int ir, ic
		double dist = 0.0

	for ir in xrange(tileA.rows):
		for ic in xrange(tileA.cols):
			dist += (tileA.mvsymbols[ir,ic] - tileB.mvsymbols[ir,ic]) ** 2
	dist = dist ** 0.5
	return dist

def getDist(taxA, taxB):
	"""
	Returns Kulczynski distance between two Tile objects (float). Tiles should
	have the same dimensions, unit polygon (geometry attribute), and not be null
	(at least a cell value is > 0).

	Arguments are two Tile objects.
	"""
	assert isinstance(taxA, Tile), "data.detDist function called on a {0} object.".format(type(taxA))
	assert isinstance(taxB, Tile), "data.detDist function called on a {0} object.".format(type(taxB))
	assert taxA.isNull() == False and taxB.isNull() == False, "data.detDist function called on a null Tile."
	assert taxA.rows == taxB.rows and taxA.cols == taxB.cols, "data.detDist function called on Tile objects of different dimensions."
	assert taxA.geometry == taxB.geometry, "data.detDist function called on Tile objects of different geometry: `{0}` and `{1}`.".format(taxA.geometry, taxB.geometry)
	return getDistance(taxA, taxB)

cdef class Tile:
	"""
	Tile objects are geographic distribution lattices. Cells values are floats
	between 0.0 and 1.0 (0.0 = absence, 1.0 = presence).

	Read-only attributes:

	- rows (int): Number of rows in the lattice.

	- cols (int): Number of cols in the lattice.

	- geometry (str, "square" or "hexagon"): geometric shape of cells making
	up the lattice.

	- name (str): a string identifier for the object, usually the taxa name.

	Instantiation arguments (__cinit__ method):

	- ingrid (list of list of ints): A two-dimensional list representing the
	geographic distribution (0 as absence, 1 as presence). If argument `template`
	is not used, `ingrid` is mandatory.

	- cellType (str, "square" or "hexagon"): Sets the basic polygon unit of the
	lattice. If If argument `template` is not used, `cellType` is mandatory.

	- template (Tile obj): Another Tile object that will serve as template for the
	new instance (dimensions and distributional data will be identical). If
	arguments `ingrid` and `celltype` are not provided, using `template` is
	mandatory.

	- name (str, optional): a string identifier for the object, usually the taxa name.
	"""
	#__name__ = "Tile"
	#cdef:
	#	double[:,::1] mvsymbols
	#	long[:,:,:,::1] mvneighref
	#	double value
	#	int ir, ic, ine, ind, ind2, neighsNum
	#	readonly int rows, cols
	#	readonly str geometry, name

	def __cinit__(self, list ingrid = None, str cellType = "square", template = None, str name = "Nameless_Ghoul", *args, **kwargs):

		if isinstance(cellType, str):
			if cellType == "square" or cellType == "hexagon":
				self.geometry = cellType
			else:
				raise ValueError("Valid arguments for \"cellType\" option are \"square\" or \"hexagon\".")

		if isinstance(name, str):
			self.name = name

		if isinstance(ingrid,list):
			if len(ingrid) >= 1 and len(ingrid[0]) >= 1 and (isinstance(ingrid[0][0],int) or isinstance(ingrid[0][0],float)):
				self.rows = len(ingrid)
				self.cols = len(ingrid[0])
				if self.geometry == "square":
					self.neighsNum = 4
				elif self.geometry == "hexagon":
					self.neighsNum = 6

				# Get memory view
				self.mvsymbols = np.empty((self.rows,self.cols), dtype=float, order = 'C')
				self.mvneighref = np.empty((self.rows,self.cols,self.neighsNum,2), dtype=long, order = 'C')
				self.mvneighref[...] = -1

				for ir in xrange(self.rows):
					for ic in xrange(self.cols):
						self.mvsymbols[ir,ic] = <double> ingrid[ir][ic]
						if self.neighsNum == 4:
							#neighs = [-1 for x in xrange(self.neighsNum)]
							ine = 0
							if ir > 0:
								self.mvneighref[ir,ic,ine,0] = ir-1
								self.mvneighref[ir,ic,ine,1] = ic
								ine += 1
							if ir < (self.rows - 1):
								self.mvneighref[ir,ic,ine,0] = ir+1
								self.mvneighref[ir,ic,ine,1] = ic
								ine += 1
							if ic > 0:
								self.mvneighref[ir,ic,ine,0] = ir
								self.mvneighref[ir,ic,ine,1] = ic-1
								ine += 1
							if ic < (self.cols - 1):
								self.mvneighref[ir,ic,ine,0] = ir
								self.mvneighref[ir,ic,ine,1] = ic+1
								ine += 1
						elif self.neighsNum == 6:
							ine = 0
							if ir > 0:
								self.mvneighref[ir,ic,ine,0] = ir-1
								self.mvneighref[ir,ic,ine,1] = ic
								ine += 1
								if ic > 0 and ir % 2 == 0: # even row index
									self.mvneighref[ir,ic,ine,0] = ir-1
									self.mvneighref[ir,ic,ine,1] = ic-1
									ine += 1
								if ic < (self.cols - 1) and ir % 2 == 1: # odd row index
									self.mvneighref[ir,ic,ine,0] = ir-1
									self.mvneighref[ir,ic,ine,1] = ic+1
									ine += 1
							if ir < (self.rows - 1):
								self.mvneighref[ir,ic,ine,0] = ir+1
								self.mvneighref[ir,ic,ine,1] = ic
								ine += 1
								if ic > 0 and ir % 2 == 0: # even row index
									self.mvneighref[ir,ic,ine,0] = ir+1
									self.mvneighref[ir,ic,ine,1] = ic-1
									ine += 1
								if ic < (self.cols - 1) and ir % 2 == 1: # odd row index
									self.mvneighref[ir,ic,ine,0] = ir+1
									self.mvneighref[ir,ic,ine,1] = ic+1
									ine += 1
							if ic > 0:
								self.mvneighref[ir,ic,ine,0] = ir
								self.mvneighref[ir,ic,ine,1] = ic-1
								ine += 1
							if ic < (self.cols - 1):
								self.mvneighref[ir,ic,ine,0] = ir
								self.mvneighref[ir,ic,ine,1] = ic+1
								ine += 1

		elif isinstance(template,Tile) or isinstance(template,HMRF):
			### instantiate Tile from Tile or HMRF classes
			self.geometry = template.geometry
			if self.geometry == "square":
				self.neighsNum = 4
			elif self.geometry == "hexagon":
				self.neighsNum = 6
			self.name = template.name
			self.rows = template.rows
			self.cols = template.cols
			self.mvsymbols = np.empty((self.rows,self.cols), dtype=float, order = 'C')
			self.mvneighref = np.array(template.getNeighsAll(), dtype=long, order = 'C')
			for ir in xrange(self.rows):
				for ic in xrange(self.cols):
					self.mvsymbols[ir,ic] = template[ir,ic]
					#otherNeighs = template.getNeighs(ir,ic)
					#for ine in xrange(self.neighsNum):
					#	self.mvneighref[ir,ic,ine,0] = otherNeighs[ine][0]
					#	self.mvneighref[ir,ic,ine,1] = otherNeighs[ine][1]

	def __add__(self, other):
		"""
		Adds two Tile objects cell-wise. Returns another Tile instance.
		"""
		if isinstance(other, Tile) or isinstance(other, HMRF):
			if self.geometry == other.geometry and self.rows == other.rows and self.cols == other.cols:
				newTile = Tile(template = self)
				for ir in xrange(self.rows):
					for ic in xrange(self.cols):
						newValue = self[ir,ic] + other[ir,ic]
						newTile.set(ir,ic,newValue)
			else:
				raise ValueError('A Tile can only be added to another Tile of the same shape.')
		else:
			raise TypeError('A Tile can only be added to another Tile object.')
		return newTile


	def __div__(self, double number):
		"""
		Divides all values in the Tile by the same number (float or int). Returns
		anther Tile object.
		"""
		if number == 0.0:
			raise ZeroDivisionError
		else:
			newTile = Tile(template = self)
			for ir in xrange(self.rows):
				for ic in xrange(self.cols):
					value = self[ir,ic] / number
					newTile.set(ir,ic,value)
			return newTile

	def __str__(self):
		"""
		Pretty tile print method. Cell values are rounded to integers.
		"""
		buff = ""
		if self.geometry == "hexagon":
			for ir in xrange(self.rows):
				if ir % 2 != 0:
					buff += "-"
				for ic in xrange(self.cols):
					buff += "{:.0f}-".format(self.mvsymbols[ir,ic])
				if ir % 2 != 0:
					buff = buff[:-1]
				buff += "\n"
		elif self.geometry == "square":
			for ir in xrange(self.rows):
				for ic in xrange(self.cols):
					buff += "{:.0f}-".format(self.mvsymbols[ir,ic])
				buff = buff[:-1]
				buff += "\n"
		return "{0}".format(buff)

	def __getitem__(self,pos):
		"""
		Simple access to a lattice value (returns float). Arguments are row and
		column indexes.
		"""
		ir, ic = pos
		if ir < self.rows and ic < self.cols:
			return self.mvsymbols[ir,ic]
		else:
			raise IndexError

	def __iter__(self):
		"""
		Iterates through the values (floats) stored in the lattice.
		"""
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				yield self.mvsymbols[ir,ic]

	def __reduce__(self):
		return (Tile, (self.toList(), self.geometry, None, self.name))

	def enumerate(self):
		"""
		Iterates through cells of the latttices. Returns row index (int), column
		index (int), and cell value (float).
		"""
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				yield ir,ic,self.mvsymbols[ir,ic]

	def getNeighs(self, ir, ic):
		"""
		Access to lattice indexes of a cell neighbors. Arguments are the row and
		column indexes of a cell. Returns a list of tuples, each tuple contains
		a pair of indexes.
		"""
		outList = []
		for ine in xrange(self.neighsNum):
			tup = (self.mvneighref[ir,ic,ine,0],self.mvneighref[ir,ic,ine,1])
			if tup != (-1,-1):
				outList.append(tup)
		return outList

	def getNeighsAll(self):
		outList = []
		for ir in xrange(self.rows):
			outList.append([])
			for ic in xrange(self.cols):
				outList[ir].append([])
				for ine in xrange(self.neighsNum):
					outList[ir][ic].append([self.mvneighref[ir,ic,ine,0],self.mvneighref[ir,ic,ine,1]])
		return outList

	def isNull(self):
		"""
		Checks if all cells in Tile are set to 0. Returns bool.
		"""
		out = True
		outbreaker = False
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				if self.mvsymbols[ir,ic] == 1.0:
					out = False
					outbreaker = True
					break
			if outbreaker:
				break
		return out

	def set(self, int row, int col, invalue):
		"""
		Set the value of a given cell.
		"""
		self.mvsymbols[row, col] = <double> invalue

	def toBits(self):
		"""
		Transforms the Tile into a string of zeros and ones---bitvector. Do not use on float Tiles.
		"""
		stringOut = ""
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				ind2 = <int>self.mvsymbols[ir,ic]
				stringOut += "{0}".format(ind2)
		return stringOut

	def toList(self):
		"""
		Transforms the Tile into a 2-dimensional list.
		"""
		listOut = [[0.0 for x in xrange(self.cols)] for x in xrange(self.rows)]
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				listOut[ir][ic] = self.mvsymbols[ir, ic]
		return listOut


cdef class HMRF(Tile):
	"""
	Hidden Markov Random field subclass. Contains field parameters and optimization algorithms.
	"""
	#cdef:
	#	readonly double gamma
	#	readonly double means[2]
	#	readonly double staDevs[2]
	#	double[:,:,::1] mvprobabilities # holds probabilities for both symbols in tuples.
	#	double[:,::1] mvaverobs  # holds observation mean by cell

	def __cinit__(self, *args, **kwargs):
		self.mvprobabilities = np.zeros((self.rows,self.cols,2), dtype=float, order = 'C')
		self.mvaverobs = np.zeros((self.rows,self.cols), dtype=float, order = 'C')
		if "gamma" in kwargs:
			if kwargs["gamma"] > 0.0:
				self.gamma = kwargs["gamma"]
			else:
				raise ValueError("Gamma parameter should be greater than zero.")
		else:
			self.gamma = 0.2 # this should always be positive, otherwise self.ufunction() and self.indPrior() can raise unwanted exception value
		if "mean0" in kwargs:
			self.means[0] = kwargs["mean0"]
		else:
			self.means[0] = 0.0
		if "mean1" in kwargs:
			self.means[1] = kwargs["mean1"]
		else:
			self.means[1] = 1.0
		if "staDev0" in kwargs:
			self.staDevs[0] = kwargs["staDev0"]
		else:
			self.staDevs[0] = 0.13
		if "staDev1" in kwargs:
			self.staDevs[1] = kwargs["staDev1"]
		else:
			self.staDevs[1] = 0.13
		super(Tile, self).__init__(*args, **kwargs)

	def __reduce__(self):
		return (HMRF, (self.toList(), self.geometry))

	cdef double indLike(self, float obs, double statein) except? -1000.0:
		cdef double energy
		cdef int state = <int>statein
		energy = (((obs - self.means[state]) ** 2)/(2.0 * (self.staDevs[state] ** 2))) + log(self.staDevs[state])
		return energy

	cdef double indPrior(self, int indRow, int indCol) except? -1.0:
		cdef double energy = 0.0
		cdef int ine, y, x
		for ine in xrange(self.neighsNum):
			y, x = self.mvneighref[indRow, indCol, ine, 0], self.mvneighref[indRow, indCol, ine, 1]
			if y >= 0 and x >= 0:
				if self.mvsymbols[indRow, indCol] != self.mvsymbols[y, x]:
					energy += self.gamma
		return energy

	cdef double logGaussPDF(self, double mean, double staDev, double observation) except? 1000:
		cdef double prob
		prob = -log(2 * M_PI)/2 - log(staDev**2)/2 - ((observation - mean)**2)/(2 * staDev**2)
		return prob

	cdef double logGaussPMF(self, double mean, double staDev, double observation) except? 1000:
		cdef double prob
		prob = -log(2 * M_PI)/2  - ((observation - mean)**2)/(2 * staDev**2)
		return prob

	cdef double logGaussScaled(self, double mean, double staDev, double observation) except? 1000:
		cdef double prob
		prob = - ((observation - mean)**2)/(2 * staDev**2)
		return prob

	cdef double ufunction(self, int indRow, int indCol) except? -1.0:
		cdef double u0 = 0.0, u1 = 0.0, uu = 0.0
		cdef int ine, y, x
		for ine in xrange(self.neighsNum):
			y, x = self.mvneighref[indRow, indCol, ine, 0], self.mvneighref[indRow, indCol, ine, 1]
			if y >= 0 and x >= 0:
				if self.mvsymbols[y,x] == 0.0:
					u0 += self.gamma
				elif self.mvsymbols[y,x] == 1.0:
					u1 += self.gamma
				if self.mvsymbols[y,x] == self.mvsymbols[indRow, indCol]:
					uu += self.gamma
		return log(exp(uu) / (exp(u0) + exp(u1)))


	cdef double bhadist(self):
		cdef double dist
		dist = 0.25 * log(0.25 * ((self.staDevs[0]**2 / self.staDevs[1]**2) + (self.staDevs[1]**2 / self.staDevs[0]**2) + 2.0)) + 0.25 * ((self.means[0] - self.means[1])**2 / (self.staDevs[0]**2 + self.staDevs[1]**2))
		if dist > 300.0:
			dist = 300.0
		return dist

	def setMeans(self, newMean0, newMean1):
		self.means[0] = newMean0
		self.means[1] = newMean1

	def setStaDevs(self, newSD0, newSD1):
		self.staDevs[0] = newSD0
		self.staDevs[1] = newSD1

	def setGamma(self, newGamma):
		self.gamma = newGamma

	def nullMe(self):
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				self.mvsymbols[ir, ic] = 0

	def getPostProbs(self):
		"""
		Returns a two-dimensional list with the posterior probability of symbol 1
		for each node in the lattice.
		"""
		fakeList = [[float for x in xrange(self.cols)] for y in xrange(self.rows)]
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				fakeList[ir][ic] = self.mvprobabilities[ir,ic,1]
		return fakeList

	def pseudoLike(self, obserIn, distr_form = 'raw'):
		"""
		Estimates Besag's pseudolikelihood. Requires a list of Tile objects (observations).
		"""
		if isinstance(obserIn, list) and isinstance(obserIn[0], Tile):
			observations = obserIn
		elif isinstance(obserIn, Tile):
			observations = [obserIn]
		else:
			raise TypeError('Argument to emea method should be either a list or a Tile object.')

		cdef double weight, pseudolikelihood, logLike, logPrior, logJoint
		cdef double mean, divisor
		cdef int inx = 0, iny = 0, ino = 0, fieldState, form

		pseudolikelihood = 0.0

		if distr_form == 'raw':
			form = 0
		elif distr_form == 'scaled':
			form = 1
		elif distr_form == 'pmf':
			form = 2

		#################################################
		# Estimate Bhattacharyya weight
		# This transformation of the distance ensures that the minimal weigth add
		# to the final pseudolikelihood value will be 0.0, and the maximal won't
		# enforce a math error of log function
		weight = 1.0 - (1.0 / exp(self.bhadist()))
		#################################################

		# Get mean observed value
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				mean = 0.0
				for ino in xrange(len(observations)):
					mean += <double>observations[ino][ir,ic]
				divisor = <double>len(observations)
				mean = mean / divisor
				self.mvaverobs[ir,ic] = mean

		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				fieldState = <int>self.mvsymbols[ir,ic]
				if form == 0:
					logLike = self.logGaussPDF(self.means[fieldState], self.staDevs[fieldState], self.mvaverobs[ir,ic])
				elif form == 1:
					logLike = self.logGaussScaled(self.means[fieldState], self.staDevs[fieldState], self.mvaverobs[ir,ic])
				elif form == 2:
					logLike = self.logGaussPMF(self.means[fieldState], self.staDevs[fieldState], self.mvaverobs[ir,ic])
				logPrior = self.ufunction(ir, ic)
				##########################################
				# Add the Bhattacharyya weight only once?
				##########################################
				logJoint = logLike + logPrior
				pseudolikelihood += logJoint
		pseudolikelihood += log(weight)

		return pseudolikelihood



	def emea(self, obserIn, debbug = False):
		"""
		Expectation-maximization algorithm based on Gibbs energy approximations.
		"""
		if isinstance(obserIn, list) and isinstance(obserIn[0], Tile):
			observations = obserIn
		elif isinstance(obserIn, Tile):
			observations = [obserIn]
		else:
			raise TypeError('Argument to emea method should be either a list or a Tile object.')

		# Get mean observed value
		cdef double mean = 0.0, divisor, num0, num1, den0, den1
		cdef int inx = 0, iny = 0, ino = 0
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				mean = 0.0
				for ino in xrange(len(observations)):
					mean += <double>observations[ino][ir,ic]
				divisor = <double>len(observations)
				mean = mean / divisor
				self.mvaverobs[ir,ic] = mean

		# 1. Set initial parameters

		#self.staDevs[0] = 0.1 * random.random()
		#self.staDevs[1] = 0.3 + (0.1 * random.random())
		#self.means[0] = 0.2 * random.random()
		#self.means[1] = 1 - (0.2 * random.random())

		# 2 and 3. Calculate likelihood energy and estimate labels
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				prevState = self.mvsymbols[ir,ic]
				currE = self.indLike(self.mvaverobs[ir,ic], self.mvsymbols[ir,ic]) + self.indPrior(ir, ic)
				if prevState == 1:
					self.mvsymbols[ir, ic] = 0
				elif prevState == 0:
					self.mvsymbols[ir, ic] = 1
				newE = self.indLike(self.mvaverobs[ir,ic], self.mvsymbols[ir,ic]) + self.indPrior(ir, ic)
				if currE < newE:
					self.mvsymbols[ir, ic] = prevState


		# 4. Estimate posterior distribution of all labels for every site
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				P0 = exp(self.logGaussPMF(self.means[0] , self.staDevs[0], self.mvaverobs[ir,ic]) + self.indPrior(ir, ic))
				P1 = exp(self.logGaussPMF(self.means[1] , self.staDevs[1], self.mvaverobs[ir,ic]) + self.indPrior(ir, ic))
				#P0 = self.staDevs[0] * exp(-(self.indLike(self.mvaverobs[ir,ic], 0) + self.indPrior(ir, ic)))
				#P1 = self.staDevs[1] * exp(-(self.indLike(self.mvaverobs[ir,ic], 1) + self.indPrior(ir, ic)))
				self.mvprobabilities[ir,ic,0] = P0 / (P0 + P1)
				self.mvprobabilities[ir,ic,1] = P1 / (P0 + P1)

		# 5. Update parameters
		(num0,num1,den0,den1) = (0.0,0.0,0.0,0.0)
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				num0 += self.mvprobabilities[ir,ic,0] * self.mvaverobs[ir,ic]
				den0 += self.mvprobabilities[ir,ic,0]
				num1 += self.mvprobabilities[ir,ic,1] * self.mvaverobs[ir,ic]
				den1 += self.mvprobabilities[ir,ic,1]
		self.means[0] = num0 / den0
		self.means[1] = num1 / den1

		if debbug:
			print "Mean 0:",self.means[0]
			print "Mean 1:",self.means[1],"\n"

		(num0,num1,den0,den1) = (0.0,0.0,0.0,0.0)
		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				num0 += self.mvprobabilities[ir,ic,0] * (self.mvaverobs[ir,ic] ** 2)
				den0 += self.mvprobabilities[ir,ic,0]
				num1 += self.mvprobabilities[ir,ic,1] * ((self.mvaverobs[ir,ic] - 1) ** 2)
				den1 += self.mvprobabilities[ir,ic,1]
		self.staDevs[0] = (num0 / den0) ** 0.5
		self.staDevs[1] = (num1 / den1) ** 0.5
		if self.staDevs[0] < 0.000001:
			self.staDevs[0] = 0.000001
		if self.staDevs[1] < 0.000001:
			self.staDevs[1] = 0.000001

		if debbug:
			print "Standard deviation 0:",self.staDevs[0]
			print "Standard deviation 1:",self.staDevs[1],"\n"


		# 6. Estimate labels again!

		for ir in xrange(self.rows):
			for ic in xrange(self.cols):
				prevState = self.mvsymbols[ir,ic]
				currE = self.indLike(self.mvaverobs[ir,ic], self.mvsymbols[ir,ic]) + self.indPrior(ir, ic)
				if prevState == 1:
					self.mvsymbols[ir, ic] = 0
				elif prevState == 0:
					self.mvsymbols[ir, ic] = 1
				newE = self.indLike(self.mvaverobs[ir,ic], self.mvsymbols[ir,ic]) + self.indPrior(ir, ic)
				if currE < newE:
					self.mvsymbols[ir, ic] = prevState
