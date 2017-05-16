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

class Ensemble(object):

	def __init__(self):
		self.aic = float()
		self.plic = float()
		self.pseudolikelihood = float()
		self.field2taxa = {}
		self.components = int()
		self.noise = []
		return None

	def __str__(self):
		buff = ""
		for infi, field in enumerate(self.field2taxa):
			buff += "\nArea {0}:\n{1}\n".format(infi,field)
			for ie, elem in enumerate(self.field2taxa[field]):
				buff += "Element {0}:\n{1}\n".format(ie, elem)
			buff += "#" * elem.cols
		return buff

	def briefStr(self):
		buff = ""
		for infi, field in enumerate(self.field2taxa):
			buff += "Area {0}:\n{1}\nElements:\n".format(infi,field)
			for ie, elem in enumerate(self.field2taxa[field]):
				buff += "{0}\n".format(elem.name)
			buff += "\n"
		return buff
