########################################################################################################################
# Name: Grain Size Distribution Tool
# Purpose: Find the average size of gravel in a stream network. Because the Q_2 value is localized, this will not work
# outside the Columbia River Basin, though modifying findQ_2() for a particular region should be relatively simple.
#
# Author: Braden Anderson
# Created: 31 March 2017
# Last Update: 18 May 2017
########################################################################################################################

import arcpy
import os
from GrainSizeReach import Reach
from math import sqrt
from timeit import default_timer as timer


def main(dem,
         flowAccumulation,
         streamNetwork,
         precipMap,
         clippingRegion,
         outputFolder,
         nValue,
         t_cValue,
         regionNumber,
         testing):
    """
    Our main function
    :param dem: The path to a DEM file
    :param streamNetwork: The path to a .shp file that contains our stream network
    :param precipMap: The path to a .shp file that contains polygons that have precipitation data
    :param huc10: The region that our stream network will be clipped to
    :param outputFolder: Where we want to put our output
    :param nValue: What value we use for our Manning coefficient. Important for our equation
    :param t_cValue: What value we use for our Shields stress coefficient. Important for our equation
    :param regionNumber: What region we use to calculate our Q_2 value
    :return: None
    """
    # TODO Implement projection checks. Look into sr.abbreviation, sr.alias,
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")  # We'll be using a bunch of spatial analysis tools

    if testing:
        arcpy.AddMessage("Currently in TestMode")

    """Creates the temporary data folder, where we'll put all our intermediate results"""
    if not os.path.exists(outputFolder+"\\temporaryData"):
        os.makedirs(outputFolder+"\\temporaryData")
    tempData = outputFolder + "\\temporaryData"

    streamSR = arcpy.Describe(streamNetwork).spatialReference
    demSR = arcpy.Describe(dem).spatialReference
    precipSR = arcpy.Describe(precipMap).spatialReference
    if streamSR.PCSName != demSR.PCSName != precipSR.PCSName != precipSR.PCSName:
        arcpy.AddError("DEM AND STREAM NETWORK USE DIFFERENT PROJECTIONS")

    """Creates our output folder, where we'll put our final results"""
    if not os.path.exists(outputFolder+"\GrainSize"):
        os.makedirs(outputFolder+"\GrainSize")
    outputDataPath = outputFolder+"\GrainSize"

    """Clips our stream network to a HUC10 region"""
    if clippingRegion != None:
        clippedStreamNetwork = tempData + "\clippedStreamNetwork.shp"
        arcpy.AddMessage("Clipping stream network...")
        arcpy.Clip_analysis(streamNetwork, clippingRegion, clippedStreamNetwork)
    else:
        clippedStreamNetwork = streamNetwork

    """Makes the reaches"""
    reachArray = makeReaches(testing, dem, flowAccumulation, clippedStreamNetwork, precipMap, regionNumber, tempData,
                             nValue, t_cValue)

    """Writes our output to a folder"""
    writeOutput(reachArray, outputDataPath, arcpy.Describe(streamNetwork).spatialReference)

    """Writes data to a text file. Delete in final build"""
    writeResults(reachArray, testing, outputDataPath)


def makeReaches(testing, dem, flowAccumulation, streamNetwork, precipMap, regionNumber, tempData, nValue, t_cValue):
    """
    Goes through every reach in the stream network, calculates its width and Q_2 value, and stores that data in a
    Reach object, which is then placed in an array
    :param testing: Bool, tells ushether or not we want to only do a small number of reaches for testing purposes.
    :param streamNetwork: The path to a .shp file that contains our stream network
    :param precipMap: The path to a .shp file that contains polygons that have precipitation data
    :param regionNumber: What region we use to calculate our Q_2 value
    :param tempData: Where we're going to put our temp data, because ArcGIS can't give us an easy way to get data from
    a raster at a point
    :param nValue: Our Manning coefficient
    :param t_cValue: Our Shields critical value
    :return: An array of reaches, with a calculated Grain size for each
    """

    reaches = []
    numReaches = int(arcpy.GetCount_management(streamNetwork).getOutput(0))
    numReachesString = str(numReaches)
    arcpy.AddMessage("Reaches to calculate: " + numReachesString)


    if flowAccumulation == None:
        arcpy.AddMessage("Calculating Drainage Area...")
        filledDEM = arcpy.sa.Fill(dem)
        flowDirection = arcpy.sa.FlowDirection(filledDEM)
        flowAccumulation = arcpy.sa.FlowAccumulation(flowDirection)
    cellSizeX = arcpy.GetRasterProperties_management(flowAccumulation, "CELLSIZEX")
    cellSizeY = arcpy.GetRasterProperties_management(flowAccumulation, "CELLSIZEY")
    cellSize = float(cellSizeX.getOutput(0)) * float(cellSizeY.getOutput(0))
    arcpy.SetProgressor("step", "Creating Reach 1 out of " + numReachesString, 0, numReaches, 1)

    arcpy.AddMessage("Creating Reach Array...")
    polylineCursor = arcpy.da.SearchCursor(streamNetwork, ['SHAPE@'])

    """If testing, only go through the loop once. Otherwise, go through every reach"""
    if testing:
        numTests = 500
        slopeTime = 0.0
        precipTime = 0.0
        flowAccTime = 0.0
        variableTime = 0.0
        wrapUpTime = 0.0
        start = timer()
        for i in range(numTests):
            arcpy.AddMessage("Creating Reach " + str(i + 1) + " out of " + str(numTests))
            row = polylineCursor.next()

            arcpy.AddMessage("Calculating Slope...")
            tempStart = timer()
            lastPointElevation = findElevationAtPoint(dem, row[0].lastPoint, tempData)
            firstPointElevation = findElevationAtPoint(dem, row[0].firstPoint, tempData)
            tempEnd = timer()
            slopeTime += (tempEnd - tempStart)
            arcpy.AddMessage("Time to calculate slope: " + str(tempEnd - tempStart) + " seconds")

            arcpy.AddMessage("Calculating Precipitation...")
            tempStart = timer()
            precip = findPrecipitation(precipMap, tempData, row[0].lastPoint)
            tempEnd = timer()
            precipTime += (tempEnd - tempStart)
            arcpy.AddMessage("Time to calculate precipitation: " + str(tempEnd - tempStart) + " seconds")

            arcpy.AddMessage("Calculating Flow Accumulation...")
            tempStart = timer()
            flowAccAtPoint = findFlowAccumulation(flowAccumulation, tempData, cellSize)
            tempEnd = timer()
            flowAccTime += (tempEnd - tempStart)
            arcpy.AddMessage("Time to calculate flow accumulation: " + str(tempEnd - tempStart) + " seconds")

            arcpy.AddMessage("Finding Variables...")
            tempStart = timer()
            slope = findSlope(row, firstPointElevation, lastPointElevation)
            width = findWidth(flowAccAtPoint, precip)
            q_2 = findQ_2(flowAccAtPoint, firstPointElevation, precip, regionNumber, tempData)
            tempEnd = timer()
            variableTime += (tempEnd - tempStart)
            arcpy.AddMessage("Time to calculate variables: " + str(tempEnd - tempStart) + " seconds")

            tempStart = timer()
            reach = Reach(width, q_2, slope, row[0])
            reach.calculateGrainSize(nValue, t_cValue)

            reaches.append(reach)
            arcpy.AddMessage("Reach " + str(i + 1) + " complete.")
            tempEnd = timer()
            wrapUpTime += (tempEnd - tempStart)
        end = timer()
        totalTime = end - start

        arcpy.AddMessage("Average time spent calculating slope: " + str(slopeTime / numTests) + " seconds")
        arcpy.AddMessage("Average time spent calculating precipitation: " + str(precipTime / numTests) + " seconds")
        arcpy.AddMessage("Average time spent calculating flow accumulation " + str(flowAccTime / numTests) + " seconds")
        arcpy.AddMessage("Average time spent calculating variables: " + str(variableTime / numTests) + " seconds")
        arcpy.AddMessage("Average time spent putting it together:" + str(wrapUpTime / numTests) + " seconds")
        arcpy.AddMessage("Average time per reach: " + str(totalTime / numTests) + " seconds")
    else:
        i = 0
        for row in polylineCursor:
            lastPointElevation = findElevationAtPoint(dem, row[0].lastPoint, tempData)
            firstPointElevation = findElevationAtPoint(dem, row[0].firstPoint, tempData)
            precip = findPrecipitation(precipMap, tempData, row[0].lastPoint)
            flowAccAtPoint = findFlowAccumulation(flowAccumulation, tempData, cellSize)

            slope = findSlope(row, firstPointElevation, lastPointElevation)
            width = findWidth(flowAccAtPoint, precip)
            q_2 = findQ_2(flowAccAtPoint, firstPointElevation, precip, regionNumber, tempData)

            reach = Reach(width, q_2, slope, row[0])
            reach.calculateGrainSize(nValue, t_cValue)

            reaches.append(reach)

            i += 1
            arcpy.SetProgressorLabel("Creating Reach " + str(i) + " out of " + numReachesString)
            arcpy.SetProgressorPosition()

    del row
    del polylineCursor

    arcpy.AddMessage("Reach Array Created.")

    return reaches


def getMinJanTemp(tempData):
    minJanTempMap = "C:\Users\A02150284\Documents\GIS Data\JanMinTemp\PRISM_tmin_30yr_normal_800mM2_01_asc.asc"
    pointLayer = tempData + "\pointJanTemp"
    arcpy.sa.ExtractValuesToPoints(tempData + "\point.shp", minJanTempMap, pointLayer)
    searchCursor = arcpy.da.SearchCursor(pointLayer + ".shp", "RASTERVALU")
    row = searchCursor.next()
    minJanTemp = row[0]
    del searchCursor
    del row
    return minJanTemp


def findQ_2(flowAccAtPoint, elevation, precip, regionNumber, tempData):
    """
    Returns the value of a two year flood event
    :param flowAccAtPoint: A float with flow accumulation at a point
    :param precip: A float with the precipitation at a point
    :param regionNumber: What region of Washington we're in. See the link in the comment below
    :return: Q_2 (float)
    """
    """These equations are based on the USGS database. To find your region, go to the following website:
    https://pubs.usgs.gov/fs/fs-016-01/ """
    if regionNumber == 1:
        q_2 = 0.35 * (flowAccAtPoint**0.923) * (precip ** 1.24)
    elif regionNumber == 2:
        q_2 = 0.09 * (flowAccAtPoint**0.877) * (precip ** 1.51)
    elif regionNumber == 3:
        q_2 = 0.817 * (flowAccAtPoint**0.877) * (precip ** 1.02)
    elif regionNumber == 4:
        q_2 = 0.025 * (flowAccAtPoint**0.880) * (precip ** 1.70)
    elif regionNumber == 5:
        q_2 = 14.7 * (flowAccAtPoint**0.815)
    elif regionNumber == 6:
        q_2 = 2.24 * (flowAccAtPoint**0.719) * (precip ** 0.833)
    elif regionNumber == 7:
        q_2 = 8.77 * (flowAccAtPoint**0.629)
    elif regionNumber == 8:
        q_2 = 12.0 * (flowAccAtPoint**0.761)
    elif regionNumber == 9:
        q_2 = 0.803 * (flowAccAtPoint**0.672) * (precip ** 1.16)
    elif regionNumber == 12:
        q_2 = 0.508 * (flowAccAtPoint ** 0.901) * ((elevation / 1000)**0.132) * (precip ** 0.926)
    elif regionNumber == 13:
        q_2 = 12.6 * (flowAccAtPoint ** 0.879) * ((elevation / 1000) ** -0.161)
    elif regionNumber == 14:
        q_2 = 9.49 * (flowAccAtPoint ** 0.903) * ((elevation / 1000)**0.055)
    elif regionNumber == 15:
        q_2 = 9.49 * (flowAccAtPoint ** 0.903) * ((elevation / 1000)**0.055)
    elif regionNumber == 16:
        q_2 = 0.000141 * (flowAccAtPoint ** 0.904) * (precip ** 3.25)
    elif regionNumber == 100:
        minJanTemp = getMinJanTemp(tempData)
        q_2 = .00013 * (flowAccAtPoint**0.8) * (precip ** 1.24) * ((minJanTemp + 273) ** 2.53)
    else:
        arcpy.AddError("Incorrect Q_2 value entered")

    q_2 /= 35.3147  # converts from cubic feet to cubic meters

    return q_2


def findWidth(flowAccAtPoint, precip):
    """
    Estimates the width of a reach, based on its drainage area and precipitation levels
    :param flowAccAtPoint: A float with flow accumulation at a point
    :param precip: A float with the precipitation at a point
    :return: Estimated width
    """
    width = 0.177 * (flowAccAtPoint ** 0.397) * (precip ** 0.453)  # This is the equation we're using to estimate width
    if width < .3:  # establishes a minimum width value
        width = .3
    return width


def findSlope(polyline, firstPointElevation, secondPointElevation):
    """
    Finds the average slope of the reach in question, given two elevations
    :param polyline: Used to find the length
    :param firstPointElevation: A float that has elevation data
    :param secondPointElevation: Another float that has elevation data
    :return: slope
    """
    length = polyline[0].length
    elevationDifference = abs(firstPointElevation - secondPointElevation)
    return elevationDifference/length


def findFlowAccumulation(flowAccumulation, tempData, cellSize):
    """
    Finds the flow accumulation at the point defined in the findPrecipitation function
    :param flowAccumulation: A raster containing flow accumulation data
    :param tempData: Where we can dump all our intermediary data points
    :param cellSize: The area of each cell
    :return: Flow accumulation at a point
    """
    arcpy.Buffer_analysis(tempData + "\point.shp", tempData + "\pointBuffer.shp", "20 Meters")
    arcpy.PolygonToRaster_conversion(tempData + "\pointBuffer.shp", "FID", tempData + "\pointBufferRaster.tif",
                                     cellsize=sqrt(cellSize))
    maxFlow = arcpy.sa.ZonalStatistics(tempData + "\pointBufferRaster.tif", "Value", flowAccumulation, "MAXIMUM")
    arcpy.sa.ExtractValuesToPoints(tempData + "\point.shp", maxFlow, tempData + "\\flowPoint")

    searchCursor = arcpy.da.SearchCursor(tempData + "\\flowPoint.shp", "RASTERVALU")
    row = searchCursor.next()
    flowAccAtPoint = row[0]
    del row
    del searchCursor
    flowAccAtPoint *= cellSize  # gives us the total area of flow accumulation, rather than just the number of cells
    flowAccAtPoint /= 1000000  # converts from square meters to square kilometers
    if flowAccAtPoint < 0:
        flowAccAtPoint = 0

    return flowAccAtPoint


def findPrecipitation(precipMap, tempData, point):
    """
    Finds the precipitation at a point
    :param precipMap: A feature class that has precipitation data at a point
    :param tempData: Where we dump our intermediary files
    :param point: Where we want to find precipitation
    :return: Precipitation
    """
    """
    sr = arcpy.Describe(precipMap).spatialReference
    arcpy.env.workspace = tempData
    arcpy.CreateFeatureclass_management(tempData, "point.shp", "POINT", "", "DISABLED", "DISABLED", sr)
    cursor = arcpy.da.InsertCursor(tempData + "\point.shp", ["SHAPE@"])
    cursor.insertRow([point])
    del cursor
    """
    pointLayer = tempData + "\point.shp"
    arcpy.Intersect_analysis([pointLayer, precipMap], "precipPoint")
    searchCursor = arcpy.da.SearchCursor(tempData + "\precipPoint.shp", "Inches")
    row = searchCursor.next()
    precip = row[0]
    precip *= 2.54  # converts to centimeters
    del row, searchCursor
    return precip


def findElevationAtPoint(dem, point, tempData):
    """
    Finds the elevation at a certain point based on a DEM
    :param dem: Path to the DEM
    :param point: The ArcPy Point that we want to find the elevation at
    :param tempData: Where we can dump all our random data points
    :return: Elevation at a point
    """
    """
    I can't find a good way to just pull the data straight from the raster, so instead, we're having to
    create the point in a layer of its own, then create another layer that has the elevation using the Extract Value
    to Points tool, then using a search cursor to get the elevation data. It's a mess, and it's inefficient, but it
    works. If anyone can find a better way, email me at banderson1618@gmail.com
    
    Testing new feature
    """
    sr = arcpy.Describe(dem).spatialReference
    arcpy.env.workspace = tempData
    arcpy.CreateFeatureclass_management(tempData, "point.shp", "POINT", "", "DISABLED", "DISABLED", sr)
    cursor = arcpy.da.InsertCursor(tempData+"\point.shp", ["SHAPE@"])
    cursor.insertRow([point])
    del cursor
    pointLayer = tempData+"\pointElevation"
    arcpy.sa.ExtractValuesToPoints(tempData+"\point.shp", dem, pointLayer)
    searchCursor = arcpy.da.SearchCursor(pointLayer+".shp", "RASTERVALU")
    row = searchCursor.next()
    elevation = row[0]
    del searchCursor
    del row
    return elevation

    # return float(arcpy.GetCellValue_management(dem, str(point.X) + " " + str(point.Y)).getOutput(0))


def writeResults(reachArray, testing, outputData):
    """This function is meant to save the results for future study"""
    testOutput = open(outputData + "\Data(readable).txt", "w")
    i = 0
    for reach in reachArray:
        i += 1
        testOutput.write("Reach " + str(i) + ":")
        testOutput.write("\nWidth: " + str(reach.width) + " meters")
        testOutput.write("\nQ_2: " + str(reach.q_2) + " cubic meters per second")
        testOutput.write("\nSlope: " + str(reach.slope))
        testOutput.write("\nGrain Size: " + str(reach.grainSize) + "\n\n")
    testOutput.close()

    inputDataFile = open(outputData + "Data(for_reach_construction).txt", "w")
    for reach in reachArray:
        inputDataFile.write("\n" + str(reach.width))
        inputDataFile.write("\n" + str(reach.q_2))
        inputDataFile.write("\n" + str(reach.slope))
        inputDataFile.write("\n" + str(reach.grainSize))
    inputDataFile.close()


def writeOutput(reachArray, outputDataPath, sr):
    arcpy.env.workspace = outputDataPath
    outputShape = outputDataPath + "\GrainSize.shp"
    tempLayer = outputDataPath + "\GrainSize_lyr"
    outputLayer = outputDataPath + "\GrainSize.lyr"
    arcpy.CreateFeatureclass_management(outputDataPath, "GrainSize.shp", "POLYLINE", "", "DISABLED", "DISABLED")
    arcpy.AddField_management(outputShape, "GrainSize", "DOUBLE")

    insertCursor = arcpy.da.InsertCursor(outputShape, ["SHAPE@", "GrainSize"])
    for reach in reachArray:
        insertCursor.insertRow([reach.polyline, reach.grainSize])
    del insertCursor

    arcpy.DefineProjection_management(outputShape, sr)

    arcpy.MakeFeatureLayer_management(outputShape, tempLayer)
    arcpy.SaveToLayerFile_management(tempLayer, outputLayer)


if __name__ == '__main__':
    main(os.sys.argv[1],
         os.sys.argv[2],
         os.sys.argv[3])
