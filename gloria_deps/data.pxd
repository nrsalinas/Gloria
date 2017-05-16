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

cdef double getDistance(Tile tileA, Tile tileB) except -1

cdef double euclidean(Tile tileA, Tile tileB) except -1

cdef class Tile:
	cdef:
		double[:,::1] mvsymbols
		long[:,:,:,::1] mvneighref
		double value
		int ir, ic, ine, ind, ind2, neighsNum
		readonly int rows, cols
		readonly str geometry, name

cdef class HMRF(Tile):

	cdef:
		readonly double gamma
		readonly double means[2]
		readonly double staDevs[2]
		double[:,:,::1] mvprobabilities # holds probabilities for both symbols in tuples.
		double[:,::1] mvaverobs  # holds observation mean by cell
		double indLike(self, float obs, double statein) except? -1000.0
		double indPrior(self, int indRow, int indCol) except? -1.0
		double logGaussPDF(self, double mean, double staDev, double observation) except? 1000
		double logGaussPMF(self, double mean, double staDev, double observation) except? 1000
		double logGaussScaled(self, double mean, double staDev, double observation) except? 1000
		double ufunction(self, int indRow, int indCol) except? -1.0
		double bhadist(self)
