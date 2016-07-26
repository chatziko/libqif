import csv
import math
from collections import defaultdict
from utils import *

def getPointsInGrid(dataset,gridNW,gridSE):

    gridP = []
    with open(dataset) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # row is a list
            if(gridSE[0] <= float(row[4]) <= gridNW[0] and gridNW[1] <= float(row[5]) <= gridSE[1]):
                gridP.append(row)
    
    return gridP # list in a list

def getPointsInCells(gridP,gridNW,gridSE,cellLength):

    # storage
    cells = defaultdict(list)  # {} = []

    gridNWLat = gridNW[0]
    gridNWLon = gridNW[1]
    gridSELat = gridSE[0]
    gridSELon = gridSE[1]
    for position in gridP:
        lat = float(position[4])
        lon = float(position[5])
        # x and y in metres
        x = haversine((lat,gridNWLon),(lat,lon)) *1000
        y = haversine((gridSELat,lon),(lat,lon)) *1000

        gridX = int(math.floor(x/cellLength))
        gridY = int(math.floor(y/cellLength))
        gridID = cantorPair(gridX,gridY)

        cells[gridID].append(position)

    return cells

def getGlobalPrior(gridsize,ptsgrid,ptscell):
    globalPrior = [0] * (gridsize*gridsize)
    pointsInCells = ptscell
    num = len(ptsgrid)

    for x in range(0,gridsize):
        for y in range(0,gridsize):
            id = cantorPair(x,y)
            if not pointsInCells[id] == []:
                globalPrior[(gridsize*x + y)] = len(pointsInCells[id]) / float(num)

    return globalPrior

def getGlobalPriorSemL(gridsize,ptsgrid,ptscell,simulation):
    globalPriorSemL = [0] * (gridsize*gridsize)
    pointsInGrid = ptsgrid
    pointsInCells = ptscell
    pointsInGridSemL = []
    for point in pointsInGrid:
        if simulation == point[3]:
            pointsInGridSemL.append(point)

    total_num_sim = len(pointsInGridSemL)

    for x in range(0,gridsize):
        for y in range(0,gridsize):
            num_in_cell_sim = 0
            id = cantorPair(x,y)
            if not pointsInCells[id] == []:
                for point in pointsInCells[id]:
                    if simulation == point[3]:
                        num_in_cell_sim += 1
                globalPriorSemL[(gridsize*x + y)] = num_in_cell_sim / float(total_num_sim)

    return globalPriorSemL


def getUserPrior(gridsize,ptsgrid,ptscell,writerFile,threshold):
    userNum = defaultdict(int) # initialize to 0
    pointsInGrid  = ptsgrid
    for i in pointsInGrid:
        userNum[int(i[0])] += 1

    out = csv.writer(open(writerFile,"w"), delimiter=',')
    pointsInCells = ptscell
    for z in userNum:
        userPrior = []
        boo = False
        for x in range(0,gridsize):
            for y in range(0,gridsize):
                id = cantorPair(x,y)
                localCount = 0
                if not pointsInCells[id] == []:                    
                    for j in pointsInCells[id]:
                        if int(j[0]) == z:
                            localCount += 1
                
                userPrior.append(localCount/float(userNum[z]))
                if localCount > threshold:
                    boo = True
        if boo:
            out.writerow(userPrior)

def getUserPriorSemL(gridsize,ptsgrid,ptscell,writerFile,threshold_sim,simulation):
    userNumSemL = defaultdict(int) # initialize to 0
    pointsInGrid  = ptsgrid
    pointsInCells = ptscell

    out = csv.writer(open(writerFile,"w"), delimiter=',')

    for point in pointsInGrid:
        if simulation == point[3]:
            userNumSemL[int(point[0])] += 1
    print len(userNumSemL)
    c = 0
    for z in userNumSemL:
        userPriorSemL = []
        boo = False
        for x in range(0,gridsize):
            for y in range(0,gridsize):
                id = cantorPair(x,y)
                localCount_sim = 0
                if not pointsInCells[id] == []:
                    for point in pointsInCells[id]:
                        if (int(point[0]) == z and simulation == point[3]):
                            localCount_sim += 1
                if localCount_sim >= threshold_sim:
                    boo = True
                userPriorSemL.append(localCount_sim/float(userNumSemL[z]))
        if boo:
            out.writerow(userPriorSemL)
        c += 1
        print c