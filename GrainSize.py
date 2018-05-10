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
from XMLBuilder import XMLBuilder
from datetime import datetime
import uuid


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
    :param clippingRegion: The region that our stream network will be clipped to
    :param outputFolder: Where we want to put our output
    :param nValue: What value we use for our Manning coefficient. Important for our equation
    :param t_cValue: What value we use for our Shields stress coefficient. Important for our equation
    :param regionNumber: What region we use to calculate our Q_2 value
    :return: None
    """
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
        raise Exception("DEM AND STREAM NETWORK USE DIFFERENT PROJECTIONS")

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
    writeOutput(reachArray, outputFolder, dem, streamNetwork, precipMap, nValue, t_cValue, regionNumber)


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
        numTests = 10
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
            i += 1
            try:
                if lastPointElevation < 0 or firstPointElevation < 0 or precip < 0 or flowAccAtPoint < 0:
                    raise Exception("Something wasn't found properly")
                slope = findSlope(row, firstPointElevation, lastPointElevation)
                width = findWidth(flowAccAtPoint, precip)
                q_2 = findQ_2(flowAccAtPoint, firstPointElevation, precip, regionNumber, tempData)

                reach = Reach(width, q_2, slope, row[0])
                reach.calculateGrainSize(nValue, t_cValue)

                reaches.append(reach)

                arcpy.SetProgressorLabel("Creating Reach " + str(i) + " out of " + numReachesString)
                arcpy.SetProgressorPosition()
            except Exception:
                if lastPointElevation < 0 or firstPointElevation < 0:
                    arcpy.AddWarning("Elevation was not found properly for reach " + str(i))
                elif precip < 0:
                    arcpy.AddWarning("Precip was not found properly for reach " + str(i))
                elif flowAccAtPoint < 0:
                    arcpy.AddWarning("Flow accumulation was not found properly for reach " + str(i))

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
    elif regionNumber == 26:
        q_2 = 4150 * (flowAccAtPoint ** 0.553) * ((elevation / 1000) ** -2.45)
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
    try:
        arcpy.Buffer_analysis(tempData + "\point.shp", tempData + "\pointBuffer.shp", "20 Meters")
        arcpy.PolygonToRaster_conversion(tempData + "\pointBuffer.shp", "FID", tempData + "\pointBufferRaster.tif",
                                     cellsize=sqrt(cellSize))    
        maxFlow = arcpy.sa.ZonalStatistics(tempData + "\pointBufferRaster.tif", "Value", flowAccumulation, "MAXIMUM")
        arcpy.sa.ExtractValuesToPoints(tempData + "\point.shp", maxFlow, tempData + "\\flowPoint")
    except:
        return -9999
    

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
    try:
        pointLayer = tempData + "\point.shp"
        arcpy.Intersect_analysis([pointLayer, precipMap], "precipPoint")
        searchCursor = arcpy.da.SearchCursor(tempData + "\precipPoint.shp", "Inches")
        row = searchCursor.next()
        precip = row[0]
        precip *= 2.54  # converts to centimeters
        del row, searchCursor
        return precip
    except:
        return -9999


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


def writeOutput(reachArray, outputFolder, dem, streamNetwork, precipMap, nValue, t_cValue, regionNumber):
    """
    Writes the output to a project folder
    :param reachArray: The array of streams with grain size values
    :param outputFolder: Where we want to put our stuff
    :param dem: The DEM given as input
    :param streamNetwork: The stream network given as input
    :param precipMap: The precipitation map given as input
    :param nValue: The number given for the Manning Coefficient
    :param t_cValue: The number given for the Critical Shields Stress
    :param regionNumber: The region we used to calculate Q_2
    :return:
    """
    projectFolder = makeFolder(outputFolder, "GrainSizeProject")
    writeInputs(projectFolder, dem, streamNetwork, precipMap, nValue, t_cValue, regionNumber)
    writeAnalyses(projectFolder, reachArray, arcpy.Describe(streamNetwork).spatialReference)
    writeXML(projectFolder)


def writeInputs(projectFolder, dem, streamNetwork, precipMap, nValue, t_cValue, regionNumber):
    """
    Copies over the inputs into the folder structure
    :param projectFolder: Where we put everything
    :param dem: The DEM given as input
    :param streamNetwork: The stream network given as input
    :param precipMap: The precipitation map given as input
    :param nValue: The number given for the Manning Coefficient
    :param t_cValue: The number given for the Critical Shields Stress
    :param regionNumber: The region we used to calculate Q_2
    :return:
    """
    inputsFolder = makeFolder(projectFolder, "01_Inputs")

    demFolder = makeFolder(inputsFolder, "01_DEM")
    demCopy = os.path.join(demFolder, os.path.basename(dem))
    arcpy.Copy_management(dem, demCopy)

    streamNetworkFolder = makeFolder(inputsFolder, "02_StreamNetwork")
    streamNetworkCopy = os.path.join(streamNetworkFolder, os.path.basename(streamNetwork))
    arcpy.Copy_management(streamNetwork, streamNetworkCopy)

    precipFolder = makeFolder(inputsFolder, "03_Precipitation")
    precipCopy = os.path.join(precipFolder, os.path.basename(precipMap))
    arcpy.Copy_management(precipMap, precipCopy)

    otherValuesFolder = makeFolder(inputsFolder, "04_OtherValues")
    with open(os.path.join(otherValuesFolder, "OtherValues.txt"), 'w') as textFile:
        textFile.write("n Value: " + str(nValue) + "\n")
        textFile.write("t_c Value: " + str(t_cValue) + "\n")
        textFile.write("Region Number: " + str(regionNumber) + "\n")



def writeAnalyses(projectFolder, reachArray, sr):
    """
    Writes the outputs
    :param projectFolder: Where we put our analyses
    :param reachArray: Our array of reaches with grain size data
    :param sr: The spatial reference of the output
    :return:
    """
    analysesFolder = makeFolder(projectFolder, "02_Analyses")
    outputFolder = getOutputFolder(analysesFolder)

    outputShape = outputFolder + "\GrainSize.shp"
    tempLayer = outputFolder + "\GrainSize_lyr"
    outputLayer = outputFolder + "\GrainSize.lyr"
    arcpy.CreateFeatureclass_management(outputFolder, "GrainSize.shp", "POLYLINE", "", "DISABLED", "DISABLED", sr)
    arcpy.AddField_management(outputShape, "GrainSize", "DOUBLE")

    insertCursor = arcpy.da.InsertCursor(outputShape, ["SHAPE@", "GrainSize"])
    for reach in reachArray:
        insertCursor.insertRow([reach.polyline, reach.grainSize])
    del insertCursor

    arcpy.MakeFeatureLayer_management(outputShape, tempLayer)
    arcpy.SaveToLayerFile_management(tempLayer, outputLayer)


def getOutputFolder(analysesFolder):
    """
    Gets us the first untaken Output folder number, makes it, and returns it
    :param analysesFolder: Where we're looking for output folders
    :return: String
    """
    i = 1
    outputFolder = os.path.join(analysesFolder, "Output_" + str(i))
    while os.path.exists(outputFolder):
        i += 1
        outputFolder = os.path.join(analysesFolder, "Output_" + str(i))

    os.mkdir(outputFolder)
    return outputFolder


def makeFolder(pathToLocation, newFolderName):
    """
    Makes a folder and returns the path to it
    :param pathToLocation: Where we want to put the folder
    :param newFolderName: What the folder will be called
    :return: String
    """
    newFolder = os.path.join(pathToLocation, newFolderName)
    if not os.path.exists(newFolder):
        os.mkdir(newFolder)
    return newFolder


def writeXML(projectFolder):
    """
    Creates an XML file for the project
    :param projectFolder: The root of the project folder
    :return: None
    """
    newXMLFilePath = os.path.join(projectFolder, 'project.rs.xml')
    if os.path.exists(newXMLFilePath):
        os.remove(newXMLFilePath)

    new_xml_file = XMLBuilder(newXMLFilePath,
                              "Project",
                              [('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance"),
                               ('xsi:noNamespaceSchemaLocation',
                                "https://raw.githubusercontent.com/Riverscapes/Program/master/Project/XSD/V1/Project.xsd")])


    new_xml_file.add_sub_element(new_xml_file.root, "ProjectType", "GrainSize")
    addMetaData(new_xml_file)
    addRealizations(new_xml_file, projectFolder)

    new_xml_file.write()


def addMetaData(new_xml_file):
    """
    Writes the metadata component of the XML file
    :param new_xml_file: The XML builder that we're using
    :return:
    """
    meta_data_element = new_xml_file.add_sub_element(new_xml_file.root, "MetaData")

    now = datetime.now()
    new_xml_file.add_sub_element(meta_data_element, "Meta", now.isoformat(), [("name", "CreatedOn")])


def addRealizations(new_xml_file, projectFolder):
    """
    Writes the Realizations part of the XML file
    :param new_xml_file: the XML Builder that we're using
    :return:
    """
    realizations_element = new_xml_file.add_sub_element(new_xml_file.root, "Realizations")
    grain_size_element = new_xml_file.add_sub_element(realizations_element,
                                                      "GrainSize",
                                                      tags=[("dateCreated", datetime.now().isoformat()),
                                                       ("guid", getUUID()),
                                                       ("productVersion", "0.0.1")])
    writeInputsXML(new_xml_file, grain_size_element, os.path.join(projectFolder, "01_Inputs"))

    analysesFolder = os.path.join(projectFolder, "02_Analyses")
    analyses_element = new_xml_file.add_sub_element(grain_size_element, "Analyses")
    for folder in os.listdir(analysesFolder):
        writeAnalysisXML(new_xml_file, analyses_element, os.path.join(analysesFolder, folder))


def writeInputsXML(new_xml_file, grain_size_element, input_folder):
    """
    Writes the inputs for the grain size project
    :param new_xml_file: The XMLBuilder
    :param grain_size_element: The XML node that we want to write to
    :param input_folder: The folder that contains all the inputs
    :return:
    """
    inputs_element = new_xml_file.add_sub_element(grain_size_element, "Inputs")

    writeSingleInputXML(new_xml_file, inputs_element, "DEM", "01_DEM", input_folder)
    writeSingleInputXML(new_xml_file, inputs_element, "StreamNetwork", "02_StreamNetwork", input_folder)
    writeSingleInputXML(new_xml_file, inputs_element, "Precipitation", "03_Precipitation", input_folder)
    writeSingleInputXML(new_xml_file, inputs_element, "OtherValues", "04_OtherValues", input_folder)


def writeSingleInputXML(new_xml_file, inputs_element, nameOfInput, folderName, input_folder):
    folderPath = os.path.join(input_folder, folderName)

    items_element = new_xml_file.add_sub_element(inputs_element, nameOfInput + "s")
    for file in os.listdir(folderPath):
        if file.endswith(".shp") or file.endswith(".tif") or file.endswith("lyr") or file.endswith(".txt"):
            item_element = new_xml_file.add_sub_element(items_element, nameOfInput)
            new_xml_file.add_sub_element(item_element, "Name", file)
            new_xml_file.add_sub_element(item_element, "Path", "01_Inputs\\" + folderName + "\\" + file)


def writeAnalysisXML(new_xml_file, analyses_element, output_folder):
    """
    Writes an analysis
    :param new_xml_file: The XMLBuilder
    :param analyses_element: The XML node that we want to write to
    :param output_folder: The folder that contains the output that we want to write XML for
    :return:
    """
    analysis_element = new_xml_file.add_sub_element(analyses_element, "Analysis")
    new_xml_file.add_sub_element(analysis_element, "Name", os.path.basename(output_folder))

    model_element = new_xml_file.add_sub_element(analysis_element, "Model")
    outputs_element = new_xml_file.add_sub_element(model_element, "Outputs")
    for file in os.listdir(output_folder):
        if file.endswith("lyr"):
            outputFolderBasename = os.path.basename(output_folder)
            stream_network_element = new_xml_file.add_sub_element(outputs_element, "StreamNetwork")
            new_xml_file.add_sub_element(stream_network_element, "Name", file)
            new_xml_file.add_sub_element(stream_network_element, "Path", "02_Analyses\\" + outputFolderBasename + "\\" + file)



def getUUID():
    return str(uuid.uuid4()).upper()


if __name__ == '__main__':
    main(os.sys.argv[1],
         os.sys.argv[2],
         os.sys.argv[3])
