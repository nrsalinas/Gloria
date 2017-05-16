#!/usr/bin/env python

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


import argparse
import os
import datetime

from gloria_deps import infile
from gloria_deps import search
from gloria_deps.gui import graph

version =  "0.3"
logfile = ""
argPass = True
today = datetime.datetime.now()
outfileRootDefault = today.strftime("Gloria_output_%Y-%m-%d_%H:%M:%S")
bufferLog = "Gloria ver. {0}\nAnalysis executed on {1}\n".format(version, today)

parser = argparse.ArgumentParser(description = 'Gloria (Geographic Location--hidden markOv Random fIeld Analysis): a Python software to delimit areas of endemism through Hidden Markov Random Fields.')

parser.add_argument('-i', '--infile', required = True, dest = 'infile', metavar = '<input_file>', action = 'store', help = 'Input file in csv format. See manual for a detailed guideline.')

parser.add_argument('-o', '--outfile_root', dest = 'outfileRoot', metavar = '<outfile_root_name>', default = None, action = 'store', help = 'Outfiles name root.')

parser.add_argument('-s', '--cell_size', required = True, dest = 'cellSize', metavar = '<#>', action = 'store', type = float, help = 'Grid cell size in geographic degrees.')

parser.add_argument('-t', '--cell_type', dest = 'cellType', metavar = 'square|hexagon', action = 'store', default = 'square', type = str, help = 'Grid cell type (`square` or `hexagon`). Default = `square`.')

parser.add_argument('-x', '--long_offset', dest = 'lonOffset', metavar = '<#>', action = 'store', default = 0, help = 'Longitudinal offset of the W border of the grid in relation to the westernmost point in the infile. Should be a positive value in geographic degrees. Default = 0.')

parser.add_argument('-y', '--lat_offset', dest = 'latOffset', metavar = '<#>', action = 'store', default = 0, help = 'Latitudinal offset of the N border of the grid in relation to the northernmost point in the infile. Should be a positive value in geographic degrees. Default = 0.')

parser.add_argument('-c', '--cohesion', dest = 'cohesion', metavar = '<#>', default = 0.3, action = 'store', type = float, help = 'Clustering cohesion parameter (0.0--1.0). Default = 0.3.')

parser.add_argument('-g', '--gamma', dest = 'gamma', metavar = '<#>', default = 20, action = 'store', type = float, help = 'Potts model gamma parameter. Default = 20.')

parser.add_argument('-m', '--max_combinations', dest = 'maxCombinations', metavar = '<#>', default = None, action = 'store', type = int, help = 'Maximum number of field combinations allowed per model optimization stage. Default = None.')

parser.add_argument('-d', '--debbug', action = 'store_true', dest = 'debbug', default = False, help = 'Executes developper\'s version.')

parser.add_argument('-v', '--version', action = 'version', version = "Gloria v. {0}".format(version))

args = parser.parse_args()

if not os.path.isfile(args.infile):
	print "\nArgument error '-i'.\n\n{0} could not be opened.\n".format(args.infile)
	argPass = False

if args.cohesion < 0.0 or args.cohesion > 1.0:
	print "\Argument error '-c'.\n\n{0} is not a valid cohesion value (should be 0.0--1.0).\n".format(args.cohesion)
	argPass = False

if args.cellType != "square" and args.cellType != "hexagon":
	print "\nArgument error '-t'.\n\n{0} is not a valid cell type value (should be `square` or `hexagon`).\n".format(args.cellType)
	argPass = False

if args.cellSize < 0:
	print "\Argument error '-s'.\n\n{0} is not a valid cell size value (should be a float or integer greater than zero).\n".format(args.cellSize)
	argPass = False

if args.lonOffset < 0:
	print "\Argument error '-x'.\n\n{0} is not a valid longitudinal offset value (should be a float or integer greater than zero).\n".format(args.lonOffset)
	argPass = False

if args.latOffset < 0:
	print "\Argument error '-y'.\n\n{0} is not a valid latitudinal offset value (should be a float or integer greater than zero).\n".format(args.latOffset)
	argPass = False

if argPass:

	if args.outfileRoot:
		if os.path.isfile("{0}.log".format(args.outfileRoot)) or os.path.isfile("{0}_basegrid.geojson".format(args.outfileRoot)) or os.path.isdir("{0}_areas_geojson_files".format(args.outfileRoot)):
			args.outfileRoot = outfileRootDefault
			print "\nWarning: output file root name changed to `{0}`".format(args.outfileRoot)
	else:
		args.outfileRoot = outfileRootDefault

	logfile = "{0}.log".format(args.outfileRoot)
	indata = infile.InputData(infile = args.infile)
	totTaxa, uniqPoints = indata.getStats()
	bufferLog += "Input data\nInfile: {0}\nTotal taxa processed: {1}\nUnique taxon-point pairs: {2}\n\n".format(args.infile, totTaxa, uniqPoints)

	bufferLog += "Grid parameters\nCell size: {0} degrees\nCell shape: {1}\nLongitudinal W offset: {2}\nLatitudinal N offset: {3}\n\nAnalysis parameters\nClustering cohesion value: {4}\nPotts model gamma: {5}\nMaximum combinations by optimization cycle: {6}\n".format(args.cellSize, args.cellType, args.lonOffset, args.latOffset, args.cohesion, args.gamma, args.maxCombinations)

	myTiles = indata.getTiles(cellSize = args.cellSize, geometry = args.cellType, offsetLat = args.latOffset, offsetLon = args.lonOffset)

	if args.cellType == "square":
		bufferLog += "Effective NW corner: {0}, {1}\n\n".format(indata.originN[0], indata.originN[1])
	elif args.cellType == "hexagon":
		pass

	result = search.fieldOptim(observations = myTiles, clusCohesion = args.cohesion, maxCycleIters = args.maxCombinations, gammaParameter = args.gamma, progress2stdout = True, debbug = args.debbug)
	bufferLog += "Results\n{0} Areas of Endemism found.\nPseudolikelihood of the prefered hypothesis: {1}\nPseudolikelihood Information Criterion (PLIC) value: {2}\nAkaike Information Criterion (AIC) value: {3}\n\n".format(result.components, result.pseudolikelihood, result.plic, result.aic)

	for indaoe, aoe in enumerate(result.field2taxa):
		bufferLog += "=" * 20
		bufferLog += "\n\nArea {0}\n\n{1}\nPosterior probabilities on area membership:\n\n".format((indaoe + 1), aoe)
		theseProbs = aoe.getPostProbs()
		for row in theseProbs:
			strss = map(lambda x: "{0:.5f}".format(x), row)
			bufferLog += "--".join(strss) + "\n"
		bufferLog += "\nTaxa:\n"
		for elem in result.field2taxa[aoe]:
			bufferLog += "{0}\n".format(elem.name)
		bufferLog += "\n"

	# Output geojson file
	basegrid = graph.Grid(indata.rows, indata.cols, indata.cellSize, indata.geometry, indata.originN, indata.originS)
	os.makedirs("{0}_areas_geojson_files".format(args.outfileRoot))
	with open("{0}_basegrid.geojson".format(args.outfileRoot), "w") as bfile:
		bfile.write(basegrid.geojson())
	basegrid.res2geojson(result, "{0}_areas_geojson_files".format(args.outfileRoot))

	now = datetime.datetime.now()
	timeDiff = now - today
	bufferLog += "Execution time: {0}\n".format(timeDiff.total_seconds())

	with open(logfile, "w") as loghandle:
		loghandle.write(bufferLog)

exit()
