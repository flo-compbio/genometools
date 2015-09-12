# Copyright (c) 2015 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import os
import csv
import bisect
import gzip

def open_plain_or_gzip(fn,mode='r'):
	try:
		gzip.open(fn).next()
	except IOError:
		return open(fn,mode)
	else:
		return gzip.open(fn,mode)

def flatten(l):
	# see http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python#comment10547502_952952
	# use incomprensible list comprehension
	return [item for sublist in l for item in sublist]

def bisect_index(a, x):
	# locate the leftmost value exactly equal to x
	i = bisect.bisect_left(a, x)
	if i != len(a) and a[i] == x:
		return i
	raise ValueError

def argsort(seq):
	# see http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
	return sorted(range(len(seq)), key=seq.__getitem__)

def argmin(seq):
	return argsort(seq)[0]

def argmax(seq):
	return argsort(seq)[-1]

def read_single(fn):
	data = []
	with open_plain_or_gzip(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			data.append(l[0])
	return data

def read_all(fn,m='r'):
	data = []
	with open_plain_or_gzip(fn,m) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			data.append(l)
	return data
