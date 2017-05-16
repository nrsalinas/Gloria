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

from math import log
sin60 = 0.8660254037844386

class Polygon(object):
	"""
	Basic lattice graph object, geometrical representation of a cell in the grid.
	"""
	def __init__(self, geometry, polySize, center):
		"""
		Arguments:

		- polySize (float): Length of the horizontal axis of the polygon. Therefore,
		if the polygon is a square, it equals to the side length; if it is an
		hexagon, it equals to apothem times two. Should be given in cartographic
		degrees.

		- center (tuple of two floats): Coordinates (longitude, latitude) of the
		center of the polygon.
		"""
		self.vertices = ()
		self.center = center
		self.size = float(polySize)
		if geometry == "square":
			half = self.size / 2.0
			NW = ((self.center[0] - half), (self.center[1] + half))
			NE = ((self.center[0] + half), (self.center[1] + half))
			SE = ((self.center[0] + half), (self.center[1] - half))
			SW = ((self.center[0] - half), (self.center[1] - half))
			self.vertices = (NW, NE, SE, SW)
		elif geometry == "hexagon":
			apothem = self.size / 2
			NNW = ((self.center[0] - apothem), (self.center[1] + (0.5 * (apothem / sin60))))
			N = (self.center[0], (self.center[1] + (apothem / sin60)))
			NNE = ((self.center[0] + apothem), (self.center[1] + (0.5 * (apothem / sin60))))
			SSE = ((self.center[0] + apothem), (self.center[1] - (0.5 * (apothem / sin60))))
			S = (self.center[0], (self.center[1] - (apothem / sin60)))
			SSW = ((self.center[0] - apothem), (self.center[1] - (0.5 * (apothem / sin60))))
			self.vertices = (NNW, N, NNE, SSE, S, SSW)
		else:
			raise ValueError("`geometry` argument can only be `square` or `hexagon`.")
		return None

	def json(self):
		"""
		Outputs json string representation of the polygon.
		"""
		out = "["
		out += ','.join(map(lambda v: '{"x":%s, "y":%s}' % (v[0], v[1]), self.vertices))
		out += "]"
		return out

class Grid(object):
	"""
	Collection of Polygon objects. Can be used to draw lattices. Requires an
	InputData object to be instantiated.
	"""
	def __init__(self, rows, cols, unitSize, unitGeometry, originN, originS = None):
		self.geometry = unitGeometry
		self.originN = originN
		self.cellSize = unitSize
		self.rows = rows
		self.cols = cols

		self.myPolys = []
		if self.geometry == "square":
			thisCenter = ()
			for ir in xrange(self.rows):
				for ic in xrange(self.cols):
					thisCenter = (self.originN[0] + (self.cellSize / 2.0) + (ic * self.cellSize), (self.originN[1] - (self.cellSize / 2.0) - (ir * self.cellSize)))
					self.myPolys.append(Polygon("square", self.cellSize, thisCenter))

		if self.geometry == "hexagon":
			thisCenter = ()
			hexSide = self.cellSize / (sin60 * 2.0)
			hexDiagonal = self.cellSize / sin60
			hexSubDiagonal = hexSide + ((hexDiagonal - hexSide) / 2)

			for ir in xrange(self.rows):
				for ic in xrange(self.cols):
					if ir % 2 == 0: # col index is even
						thisCenter = ((self.originN[0] + (self.cellSize / 2.0) + (ic * self.cellSize)), (self.originN[1] - (hexSide / 2.0) - (ir * hexSubDiagonal)))
					else: # col index is odd
						thisCenter = ((self.originN[0] + self.cellSize + (ic * self.cellSize)), (self.originN[1] - (hexSide / 2.0) - (ir * hexSubDiagonal)))
					self.myPolys.append(Polygon("hexagon", self.cellSize, thisCenter))
		return None

	def geojson(self):
		buff = '{\t"type": "Feature",\n\t"geometry": { "type": "MultiPolygon",\n\t\t"coordinates": [[\n'
		for cell in self.myPolys:
			buff += '\t\t[['
			buff += '],['.join(map(lambda v: "{0}, {1}".format(v[0], v[1]), cell.vertices))
			buff += '],[{0}, {1}]]'.format(cell.vertices[0][0], cell.vertices[0][1])
			buff += ',\n'

		buff = buff.rstrip(',\n')
		buff += '\n\t\t]]\n\t\t}\t\t},\n'
		return buff

	def res2geojson(self, res, path = ""):
		pref = ""
		if path != "":
			path += "/"
		if res.components >= 1000:
			pref = "000"
		elif res.components >= 100:
			pref = "00"
		elif res.components >= 10:
				pref = "0"
		for indx, area in enumerate(res.field2taxa):
			if pref and log((indx + 1), 10) % 1 == 0:
				pref = pref[1:]
			with open("{0}area_{1}{2}.geojson".format(path, pref, (indx + 1)), "w") as areaFile:
				#buff = '{\t"type": "FeatureCollection",\n\t"features":[\n'
				buff = '\t{\t"type": "Feature",\n\t\t"geometry": { "type": "MultiPolygon",\n\t\t\t"coordinates": [[\n'
				for r in xrange(area.rows):
					for c in xrange(area.cols):
						if area[r,c] == 1.0:
							buff += '\t\t\t[['
							buff += '],['.join(map(lambda v: "{0}, {1}".format(v[0], v[1]), self.myPolys[((r * area.cols) + c)].vertices))
							buff += '],[{0}, {1}]]'.format(self.myPolys[((r * area.cols) + c)].vertices[0][0], self.myPolys[((r * area.cols) + c)].vertices[0][1])
							buff += ',\n'
							#buff += self.myPolys[((r * c) + c)].json()

				buff = buff.rstrip(',\n')
				buff += '\n\t\t\t]]\n\t\t\t},\n'
				buff += '\t\t"properties": {\n\t\t\t"Endemic_species": ["'
				buff += '", "'.join(map(lambda sp: sp.name , res.field2taxa[area]))
				buff += '"],\n\t\t\t"Area_index": "Area %s" }\n' % (indx + 1)
				buff += '\t\t}\n'
				#buff += '\t\t},\n'
				areaFile.write(buff)

		return None


	def d3me(self, mycsvfile):
		if len(self.myPolys) > 1:
			bufferito = """<!DOCTYPE html>
			<html lang="en">
				<head>
					<meta charset="utf-8">
					<title>D3 Page Template</title>
					<script type="text/javascript" src="gloria_deps/gui/d3.v3.min.js"></script>
					<style type="text/css">
					/* No style rules here yet */
					</style>
				</head>
				<body>
					<script type="text/javascript">
						var dataset;

						var wi = 960;
						var he = 1000;
						var projection = d3.geo.mercator().center([-76, 15]).scale(950);

						var path = d3.geo.path().projection(projection);

						var svg = d3.select("body")
										.append("svg")
										.attr("width", wi)
										.attr("height", he);

						d3.json("gloria_deps/gui/countries.geojson1", function(json){
																svg.selectAll("path")
																	.data(json.features)
																	.enter()
																	.append("path")
																	.attr("d", path)
																	.style("fill", "blue")
																	.style("opacity","0.3")
																	.style("stroke", "blue");
																}
							);

						d3.csv("%s", function(data){
								svg.selectAll("circle")
									.data(data)
									.enter()
									.append("circle")
									.attr("cx", function(d){return projection([d.Longitude, d.Latitude])[0]})
									.attr("cy", function(d){return projection([d.Longitude, d.Latitude])[1]})
									.attr("r", 3)
									.style("fill", "red")
									.on("click", function(d){console.log(d);});

							});
						""" % mycsvfile

			#print "{0},{1}".format(self.originN, self.originS)

			bufferito += """svg.append("circle")
							.attr("cx", function(d){return projection([%s, %s])[0]})
							.attr("cy", function(d){return projection([%s, %s])[1]})
							.attr("r", 3)
							.style("fill", "black");

						""" % (self.originN[0],self.originN[1],self.originN[0],self.originN[1])

			if self.originS:
				bufferito += """svg.append("circle")
							.attr("cx", function(d){return projection([%s, %s])[0]})
							.attr("cy", function(d){return projection([%s, %s])[1]})
							.attr("r", 3)
							.style("fill", "black");

						""" % (self.originS[0],self.originS[1],self.originS[0],self.originS[1])
			bufferito += "var poly = ["
			bufferito += ",".join(map(lambda x: x.json(), self.myPolys))
			bufferito += """];
					var polysvg = svg.selectAll("polygon")
						.data(poly)
						.enter()
						.append("polygon");

					polysvg.attr("points", function(d){
											   		return d.map(function(d){
											         			return projection([d.x, d.y]).join(",");
																	}).join(" ");
											   		}
									)
						.style("fill", "yellow")
						.style("opacity","0.3")
						.style("stroke","black");

				</script>
			</body>
		</html>"""

			with open("plot.html","w") as fhandle:
				fhandle.write(bufferito)

		else:
			raise ValueError("Grid does not have data yet.")

		return None
