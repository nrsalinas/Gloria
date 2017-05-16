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

import unittest
from .. import search
from .. import sim
from .. import data

list0 = [[1,1,1,1,1,0,0,0,0,0]] * 5 + [[0 for x in xrange(10)]] * 5
list1 = [[0 for x in xrange(10)]] * 5 + [[1,1,1,1,1,0,0,0,0,0]] * 5
list2 = [[0,0,0,0,0,1,1,1,1,1]] * 5 + [[0 for x in xrange(10)]] * 5
zeros = [[0 for x in xrange(10)] for x in xrange(10)]


class TestDataCase(unittest.TestCase):
	"""
	Testing data classes: Tile and HMRF.
	"""
	def testTileInstance(self):
		"""
		Can Tiles be instantiated?
		"""
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")
		sq00 = data.Tile(template = sq0)
		self.assertIsInstance(sq0, data.Tile, "Could not instantiate Tile class from list.")
		self.assertTrue(hex0.geometry == "hexagon", "Could not instantiate hexagonal Tile class from list.")
		self.assertIsInstance(sq00, data.Tile, "Could not instantiate Tile class from another Tile object.")
		self.assertTrue(sq0.name == "Nameless_Ghoul", "Tile name could not be set.")
		self.assertTrue((sq0.cols, sq0.rows) == (10, 10), "Square Tile dimensions are not correct.")
		self.assertTrue((hex0.cols, hex0.rows) == (10, 10), "Hexagonal Tile dimensions are not correct.")

	def testNeighborhoodSystem(self):
		"""
		Neighbors indexes are correctly assigned?
		"""
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")
		self.assertTrue([(1, 0), (0, 1)] == sq0.getNeighs(0,0), "Could not assign neighbors in square Tile")
		self.assertTrue([(8, 9), (9, 8)] == sq0.getNeighs(9,9), "Could not assign neighbors in square Tile")
		self.assertTrue([(3, 4), (5, 4), (4, 3), (4, 5)] == sq0.getNeighs(4,4), "Could not assign neighbors in square Tile")
		self.assertTrue([(4, 0), (6, 0), (5, 1)] == sq0.getNeighs(5,0), "Could not assign neighbors in square Tile")

		self.assertTrue([(1, 0), (0, 1)] == hex0.getNeighs(0,0), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(8, 9), (9, 8)] == hex0.getNeighs(9,9), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(1, 9), (1, 8), (0, 8)] == hex0.getNeighs(0,9), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(8, 0), (8, 1), (9, 1)] == hex0.getNeighs(9,0), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(3, 4), (3, 3), (5, 4), (5, 3), (4, 3), (4, 5)] == hex0.getNeighs(4,4), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(4, 0), (4, 1), (6, 0), (6, 1), (5, 1)] == hex0.getNeighs(5,0), "Could not assign neighbors in hexagonal Tile")
		self.assertTrue([(1, 5), (1, 4), (0, 4), (0, 6)] == hex0.getNeighs(0,5), "Could not assign neighbors in hexagonal Tile")

	def testString(self):
		sqst = "1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0\n"

		hexst = "1-1-1-1-1-0-0-0-0-0-\n-1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0-\n-1-1-1-1-1-0-0-0-0-0\n1-1-1-1-1-0-0-0-0-0-\n-0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0-\n-0-0-0-0-0-0-0-0-0-0\n0-0-0-0-0-0-0-0-0-0-\n-0-0-0-0-0-0-0-0-0-0\n"

		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")

		self.assertTrue(sqst == "{}".format(sq0), "__str__ method fails for square Tiles.")
		self.assertTrue(hexst == "{}".format(hex0), "__str__ method fails for haxagonal Tiles")

	def testGetterself(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")
		for ir in xrange(sq0.rows):
			for ic in xrange(sq0.cols):
				self.assertTrue(list0[ir][ic] == sq0[ir,ic],"__getitem__ method fails in square Tiles.")
		for ir in xrange(hex0.rows):
			for ic in xrange(hex0.cols):
				self.assertTrue(list0[ir][ic] == hex0[ir,ic],"__getitem__ method fails in hexagonal Tiles.")

	def testIterator(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")

		buf0, buf1 = "", ""
		for row in list0:
			for num in row:
				buf0 += "{}".format(num)
		for num in sq0:
			buf1 += "{}".format(int(num))
		self.assertTrue(buf0 == buf1, "__iter__ methods fails in square Tiles.")
		buf1 = ""
		for num in hex0:
			buf1 += "{}".format(int(num))
		self.assertTrue(buf0 == buf1, "__iter__ methods fails in hexagonal Tiles.")

	def testIsNull(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		self.assertFalse(sq0.isNull(), "isNull method fails.")
		sq0 = data.Tile(ingrid = zeros, cellType = "square")
		self.assertTrue(sq0.isNull(), "isNull method fails.")

	def testSet(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		sq0.set(0,0,0)
		self.assertTrue(0.0 == sq0[0,0], "Cannot set new values for individual cells of Tile objects.")

	def testEnumerator(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		eir = 0
		eic = 0
		for ir,ic,value in sq0.enumerate():
			self.assertTrue(eir == ir, "enumerate method cannot raise appropriate row index.")
			self.assertTrue(eic == ic, "enumerate method cannot raise appropriate col index.")
			self.assertTrue(list0[eir][eic] == int(value), "enumerate method cannot raise appropriate value.")
			if eic < (len(list0[0]) - 1):
				eic += 1
			else:
				eic = 0
				eir += 1

	def testBitter(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		hex0 = data.Tile(ingrid = list0, cellType = "hexagon")
		bit0 = reduce(lambda x,y : reduce(lambda a,b: str(a) + str(b), x) + reduce(lambda a,b: str(a) + str(b), y), list0)
		self.assertTrue(bit0 == sq0.toBits(),"toBits method cannot transform output to string of ones and zeros.")
		self.assertTrue(bit0 == hex0.toBits(),"toBits method cannot transform output to string of ones and zeros.")

	def testToList(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		newList0 = sq0.toList()
		hex1 = data.Tile(ingrid = list1, cellType = "hexagon")
		newList1 = hex1.toList()

		temp0 = [[0.0 for x in range(10)] for x in xrange(10)]
		temp1 = [[0.0 for x in range(10)] for x in xrange(10)]
		for ir in xrange(len(list0)):
			for ic in xrange(len(list0[0])):
					temp0[ir][ic] = float(list0[ir][ic])
					temp1[ir][ic] = float(list1[ir][ic])

		self.assertTrue(newList0 == temp0, "Tile object cannot be transformed into list.")
		self.assertTrue(newList1 == temp1, "Tile object cannot be transformed into list.")

	def testHMRFsubclass(self):
		sf0 = data.HMRF(ingrid = list0, cellType = "square", gamma = 1.2, mean0 = 0.01, mean1 = 0.99, staDev0 = 0.3, staDev1 = 0.35)
		hf0 = data.HMRF(ingrid = list0, cellType = "hexagon", gamma = 1.2, mean0 = 0.01, mean1 = 0.99, staDev0 = 0.3, staDev1 = 0.35)
		self.assertIsInstance(sf0, data.HMRF, "Cannot instantiate HMRF subclass.")
		self.assertIsInstance(hf0, data.HMRF, "Cannot instantiate HMRF subclass.")

	def testPseudoLikelihood(self):
		obs = []
		for i in xrange(4):
			obs.append(data.Tile(ingrid = list0, cellType = 'square'))
		obs.append(data.Tile(ingrid = list1, cellType = 'square'))
		sf0 = data.HMRF(ingrid = list1, cellType = 'square')
		psl = sf0.pseudoLike(obs)
		#print "psl=",psl
		self.assertTrue(-877.3715 == round(psl, 4), "Pseudolikelihood estimation failed for square Tile objects.")

		obs = []
		for i in xrange(4):
			obs.append(data.Tile(ingrid = list0, cellType = 'hexagon'))
		obs.append(data.Tile(ingrid = list1, cellType = 'hexagon'))
		sf0 = data.HMRF(ingrid = list1, cellType = 'hexagon')
		psl = sf0.pseudoLike(obs)
		#print "psl=",psl
		self.assertTrue(-870.4487 == round(psl, 4), "Pseudolikelihood estimation failed for hexagonal Tile objects.")

	def testGetDistance(self):
		sq0 = data.Tile(ingrid = list0, cellType = "square")
		sq1 = data.Tile(ingrid = list1, cellType = "square")
		sq2 = data.Tile(ingrid = list0, cellType = "square")
		self.assertTrue(data.getDist(sq0, sq1) == 1.0, "Kulczynski distance function fails.")
		self.assertTrue(data.getDist(sq0, sq2) == 0.0, "Kulczynski distance function fails.")

	def testEM(self):
		st0 = data.Tile(ingrid = list0, cellType = "square")
		st1 = data.Tile(ingrid = list1, cellType = "square")
		sf0 = data.HMRF(template = st0)
		sf0.emea(st1)
		self.assertTrue(sf0.toList() == st1.toList(), "EM algorithm fails in square tiles.")

		ht0 = data.Tile(ingrid = list0, cellType = "hexagon")
		ht1 = data.Tile(ingrid = list1, cellType = "hexagon")
		hf0 = data.HMRF(template = ht0)
		hf0.emea(ht1)
		self.assertTrue(hf0.toList() == ht1.toList(), "EM algorithm fails in hexagonal tiles.")

class TestSearchCase(unittest.TestCase):
	"""
	Testing search module.
	"""
	def testFieldPreSampler(self):
		st0 = data.Tile(ingrid = list0, cellType = "square")
		st1 = data.Tile(ingrid = list1, cellType = "square")
		st2 = data.Tile(ingrid = list2, cellType = "square")
		stBits = [st0.toBits(), st1.toBits(), st2.toBits()]
		obs = [data.Tile(ingrid = list0, cellType = "square") for x in xrange(5)]
		obs += [data.Tile(ingrid = list1, cellType = "square") for x in xrange(5)]
		obs += [data.Tile(ingrid = list2, cellType = "square") for x in xrange(5)]
		prefields = search.fieldPreSampler(obs, 0.2)
		for fie in prefields:
			self.assertTrue(fie.toBits() in stBits, "HMRF presampler function fails.")

	def testReduceAreas(self):
		dicto = {x:((x - 7) ** 2) for x in xrange(15)}
		suppOut = [0, 14, 1, 13, 2, 12]
		realOut = search.reduceAreas(dicto, 10, 3)
		self.assertTrue(suppOut == realOut, "Field reducer function fails (search.reduceAreas).")

	def testFieldOptim(self):
		dat = sim.getFake(0, num=3, clus=3, exun=1, inun=1, noise=2)
		res = search.fieldOptim(dat)
		self.assertTrue(res.components == 3, "Correct number of components in field optimization could not be estimated.")
		self.assertTrue(len(res.noise) == 2, "Non-clustering distributions could not be identified as such during field optimization.")
		self.assertTrue(len(res.field2taxa) == 3, "Correct number of areas could not be estimated during field optimization.")
		for fie in res.field2taxa:
			self.assertTrue(len(res.field2taxa[fie]) == 3, "Correct number of elements per area could not be estimated during field optimization.")
			for ele in res.field2taxa[fie]:
				self.assertTrue(fie.toBits() == ele.toBits(), "Distributions could not be clustered to their parent area during field optimization.")
		#self.assertTrue(0 == 1, "Testing the test module.")


if __name__ == "__main__":
	unittest.main()
