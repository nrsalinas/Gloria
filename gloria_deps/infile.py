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

import csv
import re
from math import ceil, sin, radians
import data

sin30 = 0.5
sin60 = sin(radians(60))
sin90 = 1.0

class InputData(object):
	"""
	Input data processor class. Class constructor requires a csv file (str) with
	three columns: longitude, latitude, and taxon name.
	"""
	def __init__(self,infile):
		self.points = {}
		self.minLatitude = 91.0
		self.maxLatitude = -91.0
		self.minLongitude = 181.0
		self.maxLongitude = -181.0
		self.cases = ((0,0,0), (0,-1,0), (0,-1,-1), (-1,-1,-1), (-1,0,-1), (-1,0,0)) # for 3- to 2-dimensional index conversion
		self.originN = None
		self.originS = None
		self.cellSize = None
		self.rows = None
		self.cols = None
		self.geometry = None
		self.csvfile = infile
		lineCounter = 0
		latCol = int()
		lonCol = int()

		with open(infile,'r') as fil:
			table = csv.reader(fil)
			for row in table:
				lineCounter += 1
				if lineCounter == 1 and (re.search("[^0-9\.\-]",row[1]) or re.search("[^0-9\.\-]",row[2])):
					if re.search("lon(gitude)*",row[1],flags=re.I) and re.search("lat(itude)*",row[2],flags=re.I):
						latCol, lonCol = 2, 1
						continue
					elif re.search("lon(gitude)*",row[2],flags=re.I) and re.search("lat(itude)*",row[1],flags=re.I):
						latCol, lonCol = 1, 2
						continue
					else:
						raise IOError("Input file `{0}`: column labels do not follow the required format (`Longitude`, `Latitude`).".format(infile))

				if len(row) > 3:
					raise IOError("Line {0} in `{1}` contains more than three columns.".format(lineCounter, infile))
				row[2] = re.sub("[\s\'\"]","",row[2])
				row[1] = re.sub("[\s\'\"]","",row[1])

				if len(row[latCol]) < 1:
					raise IOError("Line {0} in `{1}` do not contain latitude data.".format(lineCounter, infile))
				if re.search("[^0-9\.\-]",row[latCol]) or float(row[latCol]) < -90 or float(row[latCol]) > 90:
					raise IOError("Line {0} in `{1}` contains invalid coordinate value(s): `{2}`.".format(lineCounter, infile, row[latCol]))

				if len(row[lonCol]) < 1:
					raise IOError("Line {0} in `{1}` do not contain longitude data.".format(lineCounter, infile))
				if re.search("[^0-9\.\-]",row[lonCol]) or float(row[lonCol]) < -180 or float(row[lonCol]) > 180:
					raise IOError("Line {0} in `{1}` contains invalid coordinate value(s): `{2}`.".format(lineCounter, infile, row[lonCol]))

				lat = float(row[latCol])
				lon = float(row[lonCol])

				if len(row[0]) > 90:
					raise IOError("Are you sure {0} is a correct taxon name? (line {0} in file `{1}`)".format(row[0],lineCounter, infile))

				if row[0] in self.points:
					self.points[row[0]][(lon,lat)] = 0
				else:
					self.points[row[0]] = { (lon,lat) : 0 }

				if self.minLatitude > lat:
					self.minLatitude = lat
				if self.maxLatitude < lat:
					self.maxLatitude = lat
				if self.minLongitude > lon:
					self.minLongitude = lon
				if self.maxLongitude < lon:
					self.maxLongitude = lon

		if len(self.points) < 2:
			raise ValueError("Input file only contain distribution data from {0} species (at least two are required).".format(len(self.points)))

		return None

	def getStats(self):
		spp = len(self.points)
		uniqPoints = 0
		for sp in self.points:
			uniqPoints += len(self.points[sp])
		return (spp, uniqPoints)

	def hexcode(self, lon, lat, DAref):
		"""
		Convert geographic coordinates into 2-axial indexes of n hexagonal grid.

		Notes on three axial indexes:
		- ADA = anti-diagonal. Starts in the origin N point. It is always
			positive.
		- DA = diagonal. Reference is also origin N point, the diagonal row
			south of originN is index 0, the one north originN is 1.
			Therefore it can with negative or positive.
		- H = horizontal. Starts in the western border of the tile. It is
			always positive.

		"""
		thisRow, thisCol = 0, 0
		DAdist = (((lat - self.originS[1]) + ((lon - self.originS[0]) / (sin60 * 2.0))) * sin60)
		ADAdist = (((self.originN[1] - lat) + ((lon - self.originN[0]) / (sin60 * 2.0))) * sin60)
		DAindex = int(DAdist / (self.cellSize / 2.0)) - DAref # Correct DAindex with reference index
		ADAindex = int(ADAdist / (self.cellSize / 2.0))
		Hindex = int((lon - self.originN[0]) / (self.cellSize / 2.0))
		ada_, da_, ha_ = map(float, [ADAindex, DAindex, Hindex])
		for ca in self.cases:
			ada = ada_ + ca[0]
			da = da_ + ca[1]
			ha = ha_ + ca[2]
			if ada + da == ha and (ada - da) % 3 == 0:
				thisRow = int((ada - da) / 3)
				if thisRow % 2.0 == 0:
					thisCol = int(ha / 2)
				else:
					thisCol = int((ha - 1) / 2)
				break
		return (thisRow, thisCol)

	def getTiles(self, cellSize, geometry = "square", offsetLat = 0.0, offsetLon = 0.0):
		"""
		Create basic data structures required for the analysis from a collection
		of distributional points. Returns a list of data.Tile objects.

		Arguments:

		- geometry (str, "square" or "hexagon"): Shape of the polygonal unit of
		the lattice.

		- cellSize (int or float): Size of the cells making up the lattice. If
		the grid is made up of squares, `cellSize` will be the side of the square.
		If it is made of hexagons, it will be the length of the shorter line
		traversing the hexagon through the center (== 2 * apothem).

		"""
		self.cellSize = float(cellSize)
		self.geometry = geometry
		tileStack = []
		self.rows, self.cols = 0, 0
		offsetLat = float(offsetLat)
		offsetLon = float(offsetLon)
		correctionFactor = self.cellSize / 100

		if geometry == "square":
			self.originS = None
			self.originN = ((self.minLongitude - offsetLon), (self.maxLatitude + offsetLat))
			span = (max((self.maxLongitude - self.originN[0]), self.cellSize), max((self.originN[1] - self.minLatitude), self.cellSize))
			totCols = int(ceil(span[0] / self.cellSize))
			totRows = int(ceil(span[1] / self.cellSize))
			self.rows, self.cols = totRows, totCols
			for taxon in self.points:
				grid = [[0 for x in xrange(totCols)] for x in xrange(totRows)]
				for lon,lat in self.points[taxon]:
					apprindx = ceil(((lon - self.originN[0]) / span[0]) * totCols)
					apprindy = ceil(((self.originN[1] - lat) / span[1]) * totRows)
					if apprindx == 0:
						x = 0
					else:
						x = int(apprindx - 1)
					if apprindy == 0:
						y = 0
					else:
						y = int(apprindy - 1)
					grid[y][x] = 1
				tileStack.append(data.Tile(ingrid = grid, cellType = geometry, name = taxon))

		elif geometry == "hexagon":
			# Get handy hexagonal dimensions
			hexSide = self.cellSize / (sin60 * 2.0)
			hexDiagonal = hexSide * 2.0
			hexSubDiagonal = hexSide + ((hexDiagonal - hexSide) / 2)
			#print "hexSide: ",hexSide,", hexDiagonal: ",hexDiagonal

			self.originN = ((self.minLongitude - offsetLon - (self.cellSize / 2.0)), (self.maxLatitude + offsetLat))
			span = (max((self.maxLongitude - self.originN[0]), self.cellSize), max((self.originN[1] - self.minLatitude), hexSide))
			#print "originN: ",self.originN
			#print "span: ",span
			#span = ((self.maxLongitude - self.originN[0]), (self.originN[1] - self.minLatitude))

			totRows, totCols = 0, 0

			# Get origin south
			totRows = int(ceil(abs(span[1] - hexSide) / hexSubDiagonal)) + 1
			if totRows % 2 == 0: # even number of rows
				self.originS = ((self.originN[0] + (self.cellSize / 2)), (self.originN[1] - hexSide - ((totRows - 1) * hexSubDiagonal)))
			else: # odd number of cols
				self.originS = (self.originN[0], (self.originN[1] - hexSide - ((totRows - 1) * hexSubDiagonal)))

			totCols = int(ceil(span[0] / self.cellSize))
			#print "originS: ",self.originS

			#print "totRows: {0}, totCols {1}\n".format(totRows, totCols)
			#debugBuffer +=  "totRows: {0}, totCols {1}\n".format(totRows, totCols)
			self.rows, self.cols = totRows, totCols


			# Find DA reference index
			if totRows % 2 == 0: # even number of rows
				DAref = int((((self.originN[1] - self.originS[1]) - ((self.originS[0] - self.originN[0]) / (sin60 * 2.0)) + (hexSide / 2.0)) * sin60) / (self.cellSize / 2.0)) - 1
			else:
				DAref = int(((self.originN[1] - self.originS[1] + (hexSide / 2.0)) * sin60) / (self.cellSize / 2.0)) - 1

			cases = ((0,0,0), (0,-1,0), (0,-1,-1), (-1,-1,-1), (-1,0,-1), (-1,0,0)) # for 3- to 2-dimensional index conversion
			for taxon in self.points:
				grid = [[0 for x in xrange(totCols)] for x in xrange(totRows)]
				for lon,lat in self.points[taxon]:
					#print "lon: ",lon,", lat: ",lat
					thisRow, thisCol = self.hexcode(lon, lat, DAref)
					try:
						grid[thisRow][thisCol] = 1
					except IndexError:
						fakeLon, fakeLat = lon, lat
						if thisCol == totCols:
							fakeLon -= correctionFactor
						elif thisCol < 0:
							fakeLon += correctionFactor
						if thisRow == totRows:
							fakeLat += correctionFactor
						elif thisRow < 0:
							fakeLat -= correctionFactor
						thisRow, thisCol = self.hexcode(fakeLon, fakeLat, DAref)
						grid[thisRow][thisCol] = 1

				tileStack.append(data.Tile(ingrid = grid, cellType = geometry, name = taxon))

		else:
			raise ValueError("Valid values for argument `geometry` are `square` and `hexagon`.")

		return tileStack
