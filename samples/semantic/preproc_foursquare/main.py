from utils import *
from foursquare_preproc import *
import pickle
import os

def run(dataset,gridNW,gridSE,cellLength,gridsize,threshold,simulation,threshold_sim):

	##### Points in Grid #####
	pointsInGrid = None
	if os.path.isfile(dataset+"-pointsInGrid.p"): 
		pointsInGrid = pickle.load(open(dataset+"-pointsInGrid.p","rb"))
		print len(pointsInGrid)
	else:
		pointsInGrid = getPointsInGrid(dataset,gridNW,gridSE)
		pickle.dump(pointsInGrid, open(dataset+"-pointsInGrid.p","wb"))

	##### Points in Cells #####
	pointsInCells = None
	if os.path.isfile(dataset+"-pointsInCells.p"):
		pointsInCells = pickle.load(open(dataset+"-pointsInCells.p","rb"))
		print len(pointsInCells)
	else:		
		pointsInCells = getPointsInCells(pointsInGrid,gridNW,gridSE,cellLength)
		pickle.dump( pointsInCells, open( dataset+"-pointsInCells.p", "wb" ) )

	##### Global Prior #####

	globalPrior = None
	if os.path.isfile(dataset+"-foursquare-global.csv"):
		print "Global prior present"
	else:
		globalPrior = getGlobalPrior(gridsize,pointsInGrid,pointsInCells)	
		out = csv.writer(open(dataset+"-foursquare-global.csv","w"), delimiter=',')
		out.writerow(globalPrior)

	##### User Prior #####
	
	if os.path.isfile(dataset+"-foursquare-user.csv"):
		print "User prior present"
	else:
		writerFile = dataset+"-foursquare-user.csv"
		getUserPrior(gridsize,pointsInGrid,pointsInCells,writerFile,threshold)

	##### Global Prior for a Semantic Label #####

	globalPriorSemL = getGlobalPriorSemL(gridsize,pointsInGrid,pointsInCells,simulation)
	out = csv.writer(open(dataset+"-foursquare-global-"+simulation+".csv","w"), delimiter=',')
	out.writerow(globalPriorSemL)
	print "Global prior for simulation done"

	##### User Prior for a Semantic Label #####

	writerFile = dataset+"-foursquare-user-"+simulation+".csv"
	getUserPriorSemL(gridsize,pointsInGrid,pointsInCells,writerFile,threshold_sim,simulation)
	print "User prior for simulation done"



if __name__ == "__main__":

	"""
	Configure all parameters below before running the script

	"""

	dataset = "temp/tky.csv"
	gridNW = (35.8670,139.4788) # north-west (lat,long)
	gridSE = (35.5100,139.9100) # south-east (lat,long)
	cellLength = 400.0 # in meters
	gridsize = 100
	threshold = 20 # for users in a cell
	
	simulation = "Dessert Shop"
	threshold_sim = 1


	run(dataset,gridNW,gridSE,cellLength,gridsize,threshold,simulation,threshold_sim)

	

