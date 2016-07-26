import sys, getopt
from utils import *

def init(argv):
	inputfile = ''
   	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])

	except getopt.GetoptError:
		print "Please use the below format"
		print 'init.py -i <inputFileName> -o <outputFileName>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'init.py -i <inputFileName> -o <outputFileName>'
			print "<inputFileName> in .txt format"
			print "<outputFilename>"
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg

	txt2csv(inputfile,outputfile)

	
if __name__ == "__main__":
	init(sys.argv[1:])
	
