import csv
import os
import math
from math import radians, cos, sin, asin, sqrt
from collections import defaultdict
import pickle
AVG_EARTH_RADIUS = 6371  # in km

def txt2csv(input_file,output_file):

	txt_file = "datasets/"+input_file # input txt file (dataset)
	directory = "temp"
	if not os.path.exists(directory):
		os.makedirs(directory)
	csv_file = "temp/"+output_file # output csv filename

	# use 'with' if the program isn't going to immediately terminate
	# so you don't leave files open
	# the 'b' is necessary on Windows
	# it prevents \x1a, Ctrl-z, from ending the stream prematurely
	# and also stops Python converting to / from different line terminators
	# On other platforms, it has no effect
	in_txt = csv.reader(open(txt_file, "rb"), delimiter = '\t')
	out_csv = csv.writer(open(csv_file, 'wb'))

	out_csv.writerows(in_txt)
def haversine(point1, point2, miles=False):
	""" Calculate the great-circle distance bewteen two points on the Earth surface.

	:input: two 2-tuples, containing the latitude and longitude of each point
	in decimal degrees.

	Example: haversine((45.7597, 4.8422), (48.8567, 2.3508))

	:output: Returns the distance bewteen the two points.
	The default unit is kilometers. Miles can be returned
	if the ``miles`` parameter is set to True.

	"""
	# unpack latitude/longitude
	lat1, lng1 = point1
	lat2, lng2 = point2

	# convert all latitudes/longitudes from decimal degrees to radians
	lat1, lng1, lat2, lng2 = map(radians, (lat1, lng1, lat2, lng2))

	# calculate haversine
	lat = lat2 - lat1
	lng = lng2 - lng1
	d = sin(lat * 0.5) ** 2 + cos(lat1) * cos(lat2) * sin(lng * 0.5) ** 2
	h = 2 * AVG_EARTH_RADIUS * asin(sqrt(d))
	if miles:
		return h * 0.621371  # in miles
	else:
		return h  # in kilometers

def cantorPair(k1, k2, safe=True):
	"""
	Cantor pairing function
	http://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
	"""
	z = int(0.5 * (k1 + k2) * (k1 + k2 + 1) + k2)
	if safe and (k1, k2) != cantorDepair(z):
		raise ValueError("{} and {} cannot be paired".format(k1, k2))
	return z

def cantorDepair(z):
	"""
	Inverse of Cantor pairing function
	http://en.wikipedia.org/wiki/Pairing_function#Inverting_the_Cantor_pairing_function
	"""
	w = math.floor((math.sqrt(8 * z + 1) - 1)/2)
	t = (w**2 + w) / 2
	y = int(z - t)
	x = int(w - y)
	# assert z != pair(x, y, safe=False):
	return (x, y)

