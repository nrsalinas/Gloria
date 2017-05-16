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


from setuptools import setup, find_packages
from distutils.core import Extension

setup(name = 'Gloria',
	version = '0.3',
	author = 'Nelson R. Salinas',
	author_email = 'nrsalinas@gmail.com',
	url = 'https://github.com/nrsalinas/gloria',
	description = 'A Python program to uncover areas of endemism using Hidden Markov Random Fields.',
	license = 'GNU GPL v. 3',
	keywords = 'biogeography bioinformatics',
	#install_requires = ['pyshp'],
	test_suite = 'nose.collector',
	tests_require = ['nose'],
	ext_modules = [Extension('gloria_deps.data', ['gloria_deps/data.c'])],
	scripts = ['gloria.py'],
	packages = ['gloria_deps','gloria_deps.gui','gloria_deps.examples','gloria_deps.test']
	)
