Gloria version 0.3 May 16th, 2017

**gloria.py** is a Python program to identify areas of endemism from geographic
distribution data, using Hidden Markov Random Fields.

Input is a csv file containing localities (geographic coordinates) of
the taxa to analyze (format specification under INSTRUCTIONS).

Output is a log file and several geojson files (basegrid and one for each area of
endemism).


## REQUIREMENTS

Currently, Gloria is only supported in Linux and OS operating systems. Either
case it is required:

* a Python 2.7 interpreter.

* a C compiler.

* attention to detail.


## INSTALLING

After download, simply execute `python setup.py test` to test the source
distribution. To fully install the program type `python setup.py install`.

Installation can alternatively be done through pip: `pip install <Gloria tar file>`.


## INSTRUCTIONS

The input file should conform to the following directions:

1. Data should be organize in three columns: taxon name, longitud, and latitude.

2. The first row should contain the headers "Taxon", "Latitude", and
"Longitude".

3. Longitude and latitude should be in decimal format, with periods used as
decimal marks.

4. Datapoints belonging to the same taxa should have the same string as "Taxon".

Example:

| Taxon_name | Latitude | Longitude |
|------------|----------|-----------|
| Sp_0       | -95.67   | 1.44      |
| Sp_0       | -96.07   | 0.84      |
| Sp_1       | -85.61   | -1.68     |
| Sp_1       | -75.87   | 4.12      |

which in raw csv format should be:

    Taxon_name,Latitude,Longitude
    Sp_0,-95.67,1.44
    Sp_0,-96.07,0.84
    Sp_1,-85.61,-1.68
    Sp_1,-75.87,4.12


## WARNINGS

Poorly curated datasets usually lead to ambiguous results.

As of version 0.3, this program does not include a graphical interface.


## COPYRIGHT INFORMATION AND LICENSE

Copyright 2016-2017 Nelson R. Salinas

This file is part of Gloria.

Gloria is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Gloria is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gloria.  If not, see <http://www.gnu.org/licenses/>.


## CONTACT

Nelson R. Salinas
nrsalinas@gmail.com


## CITATION

If you use this program, please cite:

Salinas, N. R. and W. C. Wheeler. _Statistical Modeling of Distribution Patterns:
a Markov Random Field Implementation and its Application on Areas of Endemism._
In preparation.
