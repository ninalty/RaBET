import arcpy
import os
import sys
import time
import arcpy
import glob
import datetime
import shutil
from operator import itemgetter
from arcpy import env
from arcpy.sa import *
import webbrowser

class RaBETAnalysisTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "RaBET Analysis Tool"
        self.description = "This tool calculates zonal statistics for user-defined zones and produces data charts for change detection analysis."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            displayName="Input Woody Canopy Cover Rasters",
            name="rasterstring",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input",
            multiValue=True)

        param1 = arcpy.Parameter(
            displayName="Output Directory",
            name="outdir",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        param2 = arcpy.Parameter(
            displayName="Output Name",
            name="outname",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Area of Interest Input Method",
            name="method",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param3.filter.type = "ValueList"
        param3.filter.list = ["Input Shapefile", "Draw Polygon on Map Interactively"]

        param4 = arcpy.Parameter(
            displayName="Draw Polygon Interactively",
            name="Draw Polygon",
            datatype="GPFeatureRecordSetLayer",
            parameterType="Required",

            direction="Input")
        param4.enabled = "false"
        param4.value = os.path.join(sys.path[0], "tooldata", "UI.shp")
        # param3.filter.type = "ValueList"
        # param3.filter.list = []

        param5 = arcpy.Parameter(
            displayName="Area of Interest Shapefile",
            name="ROIshp",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        param5.enabled = "false"

        param6 = arcpy.Parameter(
            displayName="Area of Interest Field Name",
            name="zonenamefield",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        param6.enabled = "false"

        param6.filter.type = "ValueList"
        param6.filter.list = []

        param7 = arcpy.Parameter(
            displayName="Chart Output of Mean Woody Canopy Cover",
            name="chartmetric",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param7.value = "true"

        param8 = arcpy.Parameter(
            displayName="Output Excel Spreadsheet",
            name="excelfile",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param8.value = "true"

        param9 = arcpy.Parameter(
            displayName="Ignore NoData in Calculations",
            name="nodata",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param9.value = "true"

        params = [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        if str(parameters[3].value) == "Input Shapefile":
            parameters[4].enabled = "false"
            parameters[5].enabled = "true"
            parameters[6].enabled = "true"
            # parameters[4].value = "null"

            if parameters[5].altered:
                # field_names = []

                # parameters[6].filter.list = field_names
                featureclass = parameters[5].value
                # parameters[8].value = str(featureclass)
                field_names = sorted([f.name for f in arcpy.ListFields(featureclass)])
                # field_names.insert(0, "No Name Field")
                parameters[6].filter.list = field_names
        elif (parameters[3].value) == "Draw Polygon on Map Interactively":
            parameters[4].enabled = "true"
            parameters[5].enabled = "false"
            parameters[6].enabled = "false"
            # parameters[4].value = os.path.join(sys.path[0], "tooldata","UI.lyr")
            parameters[5].value = os.path.join(sys.path[0], "tooldata", "UI.shp")

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Version and update information
        version = "1.0"
        modified = r"4/27/2016"

        # Define input parameters
        rasterstring = parameters[0].valueAsText.replace("'", "")
        rasterlist = sorted(rasterstring.split(";"))  # Create and sort list of input rasters
        outdirmain = parameters[1].valueAsText  # Output directory - all files output here
        outname = parameters[2].valueAsText  # rootname for output files
        method = parameters[3].valueAsText  #
        featureset = parameters[4].valueAsText  #
        zoneshp = parameters[5].valueAsText  # Shapefile input
        zonename = parameters[6].valueAsText  # Field name containg unique display name

        arcpy.AddMessage(method)
        charttype = parameters[7].valueAsText  # Chart output
        excelfile = parameters[8].valueAsText  # excel sheet output
        nodata = parameters[9].valueAsText  # NoData handling in zonal statistics

        outdirname = "Analysis_" + outname
        outdir = os.path.join(outdirmain, outdirname)
        if not arcpy.Exists(outdir):
            arcpy.management.CreateFolder(outdirmain, outdirname) # TL
        else:
            arcpy.AddError("Output directory exists. Please enter a different output name.")
            sys.exit(0)

        if nodata == "true":
            datatype = "DATA"

        else:
            datatype = "NODATA"

        # Check out the Spatial Analyst extension
        arcpy.CheckOutExtension("Spatial")

        # Determine directories containg Python toolbox and reource data
        rootpath = sys.path[0]
        scratchws = os.path.join(rootpath, r"scratch\scratch.gdb")
        tooldatapath = os.path.join(rootpath, r"ToolData")
        mlradata = os.path.join(tooldatapath, r"MLRA_Data.shp")

        # Set Geoprocessing environments
        # arcpy.env.scratchWorkspace = scratchws
        arcpy.env.overwriteOutput = True
        # arcpy.env.workspace = outdir
        landsat = arcpy.sa.Raster(rasterlist[0])
        spatial_ref = arcpy.Describe(landsat).spatialReference
        arcpy.env.outputCoordinateSystem = spatial_ref

        # Output names
        outtablename = outname + ".dbf"
        outtable = os.path.join(outdir, outtablename)  # Table output name
        outshapename = outname + "_polygons"
        outshape = os.path.join(outdir, outshapename + ".shp")
        outlayer = os.path.join(outdir, outshapename + ".lyr")
        zonelayer = os.path.join("in_memory", outname + "_zonelayer")
        zonethresh = 10  # max number of graphs that can dispayed before visualization is disabled.

        # Copy input shapefile to outputdirectory and create layer in memory
        if method == "Input Shapefile":
            arcpy.management.CopyFeatures(zoneshp, outshape)
            arcpy.management.MakeFeatureLayer(outshape, zonelayer)

            if zonename != None:
                zonefield = zonename
                arcpy.AddMessage("UI fieldname: " + zonefield)
                TitleString = ""
            else:
                fieldnames = [field.name for field in arcpy.ListFields(zonelayer)]
                zonefield = fieldnames[0]  # Set the Zone Field to the first OID field
                arcpy.AddMessage("UI fieldname: " + zonefield)
                TitleString = "Polygon "

            # arcpy.arcpy.MakeFeatureLayer_management(outshp, zonelayer)


        elif method == "Draw Polygon on Map Interactively":
            arcpy.AddMessage("feature set")

            arcpy.management.CopyFeatures(featureset, outshape)
            arcpy.management.MakeFeatureLayer(outshape, zonelayer)
            arcpy.management.SaveToLayerFile(zonelayer, outlayer)

            fieldnames = [field.name for field in arcpy.ListFields(zonelayer)]
            zonefield = fieldnames[0]  # Set the Zone Field to the first OID field
            arcpy.AddMessage("UI fieldname: " + zonefield)

            TitleString = "Polygon "
        else:
            arcpy.AddMessage("Area of Interest Method Not Supported")

        # Add Zone ID field
        arcpy.management.AddField(zonelayer, "Zone_ID", "TEXT")
        expression = r"!" + zonefield + r"!"
        arcpy.management.CalculateField(zonelayer, "Zone_ID", expression, "PYTHON_9.3")

        # Get geometry of MLRA for intersection check
        firstraster = str(os.path.basename(rasterlist[0]))
        mlraID = firstraster.split("_")[2]
        arcpy.AddMessage(mlraID)
        expression = """"MLRA_ID" = '{0}'""".format(mlraID)
        cliplayertile = os.path.join("in_memory", "Clip_Layer_" + outname)
        arcpy.management.MakeFeatureLayer(mlradata, cliplayertile, expression)

        # Check for more than one MLRA in WC rasters
        MLRAList = []
        for rasters in rasterlist:
            MLRAList.append(os.path.basename(rasters).split("_")[2])
        UniqueMLRAs = set(MLRAList)
        arcpy.AddMessage(UniqueMLRAs)
        if len(UniqueMLRAs) > 1:
            arcpy.AddError(
                "WC rasters from multiple MLRAs were input. All input WC rasters must be from the same MLRA.")
            sys.exit(0)

        # Check for intersection of MLRA  boundary and area of interest polygons
        inttest = os.path.join("in_memory", "inttest")
        arcpy.analysis.Intersect([cliplayertile, zonelayer], inttest)
        polyrows = int(arcpy.management.GetCount(zonelayer).getOutput(0))
        rowcount = int(arcpy.management.GetCount(inttest).getOutput(0))
        if rowcount < polyrows:
            arcpy.AddError(
                "Error subsetting MLRA " + mlraID + "." + "\n" + "Ensure that the area of interest polygons are located inside the MLRA.")
            sys.exit(0)

        # Add RaBET_ID field
        fields = arcpy.ListFields(zonelayer)
        fieldnames = []
        for names in fields:
            fieldnames.append(str(names.name))
        if "RaBET_ID" not in fieldnames:
            arcpy.management.AddField(zonelayer, "RaBET_ID", "TEXT")

        if "RaBET_NUM" not in fieldnames:
            arcpy.management.AddField(zonelayer, "RaBET_NUM", "LONG")

        if "RaBET_NAME" not in fieldnames:
            arcpy.management.AddField(zonelayer, "RaBET_NAME", "TEXT")

        # create unique zone id with cursor
        cursor = arcpy.UpdateCursor(zonelayer)
        for id, row in enumerate(cursor):
            RaBET_ID = "Area_" + str(id + 1).zfill(2)
            RaBET_NUM = id + 1

            row.setValue("RaBET_ID", RaBET_ID)
            row.setValue("RaBET_NUM", RaBET_NUM)
            if method == "Input Shapefile" and zonename != None:

                polyname = row.getValue(zonefield)
                row.setValue("RaBET_NAME", str(polyname))

            else:
                row.setValue("RaBET_NAME", TitleString + str(RaBET_NUM))

            cursor.updateRow(row)
            del row
        del cursor

        zonelist = [zone[0] for zone in arcpy.da.SearchCursor(zonelayer, "RaBET_ID")]
        del zone
        uniquezones = sorted(set(zonelist))
        zonenum = len(uniquezones)

        # Calculate zonal statistics for each polygon for each year(raster)
        tablelist = []
        for index, zones in enumerate(uniquezones):
            expression = """"RaBET_ID" = '{0}'""".format(zones)
            arcpy.management.SelectLayerByAttribute(zonelayer, "NEW_SELECTION", expression)

            namecursor = str([row[0] for row in arcpy.da.SearchCursor(zonelayer, "RaBET_NAME")][0])

            arcpy.AddMessage("\n" + "Processing statistics for: {0}".format(str(namecursor)))

            for index, rasters in enumerate(rasterlist):
                rastername = os.path.splitext(os.path.basename(rasters))[0]
                Yearname = rastername[-4:]

                currenttable = os.path.join("in_memory", str(zones) + Yearname[-2:])

                try:

                    ZonalStatisticsAsTable(zonelayer, "RaBET_ID", rasters, currenttable, datatype, "ALL")
                    arcpy.management.AddField(currenttable, "CV", "DOUBLE")
                    # arcpy.AddField_management(currenttable, "Difference", "DOUBLE")
                    arcpy.management.AddField(currenttable, "Image", "TEXT")
                    arcpy.management.AddField(currenttable, "Year", "TEXT")
                    arcpy.management.AddField(currenttable, "RaBET_Name", "TEXT")

                    arcpy.management.DeleteField(currenttable, "ZONE-CODE")

                    cursor = arcpy.UpdateCursor(currenttable)
                    for row in cursor:
                        row.setValue('Image', rastername)
                        row.setValue('Year', Yearname)
                        row.setValue('RaBET_NAME', namecursor)

                        mean = row.getValue('MEAN')

                        std = row.getValue('STD')
                        if mean != 0:
                            CV = float(std / mean)
                            row.setValue("CV", CV)

                        cursor.updateRow(row)

                        del row
                    del cursor
                    tablelist.append(currenttable)
                    arcpy.AddMessage(Yearname)
                except:
                    arcpy.AddMessage(Yearname + " - Insufficient data")

        # Define fields to be included in merged table
        fieldMappings = arcpy.FieldMappings()
        removetables = []

        arcpy.AddMessage("\n" + "Preparing data for output." + "\n")

        for tables in tablelist:

            rowcount = str(int(arcpy.management.GetCount(tables).getOutput(0)))

            if rowcount != '1':
                removetables.append(tables)
            else:
                fieldMappings.addTable(tables)

        # Exclude tables with no data in merge
        for nulltables in removetables:
            tablelist.remove(nulltables)

        # Merge zonal statistics tables in to final output table
        allfieldslist = []
        for field in fieldMappings.fields:
            allfieldslist.append(field.name)

        if field.name not in allfieldslist:
            fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex(field.name))

        arcpy.management.Merge(tablelist, outtable, fieldMappings)

        meanvalues = []
        differencevalues = []
        titlevalues = []

        templayer = outname + "_polygons"
        arcpy.management.MakeFeatureLayer(outshape, templayer)
        arcpy.management.SaveToLayerFile(templayer, outlayer, "RELATIVE")
        arcpy.management.Delete(zonelayer)

        # Calculate difference
        for zones in uniquezones:
            count = 0
            cursor = arcpy.UpdateCursor(outtable)
            for row in cursor:
                currentzone = row.getValue("RaBET_ID")

                mean = row.getValue("MEAN")

                if currentzone == zones:
                    if count == 0:
                        difference = 0

                        initialmean = mean
                    else:
                        difference = mean - initialmean

                    # row.setValue("Difference", float(difference))

                    count += 1
                cursor.updateRow(row)
                meanvalues.append(mean)
                differencevalues.append(abs(difference))
                del row
            del cursor

        # Output Excel spreadsheet
        if excelfile == "true":
            try:
                arcpy.AddMessage("\n" + "Exporting Excel spreadsheet." + "\n")
                excelout = os.path.join(outdir, outname + ".xls")
                arcpy.conversion.TableToExcel(outtable, excelout)
            except:
                arcpy.AddError("Error exporting Excel spreadsheet.")
                pass

        # Output Chart
        if charttype == "true":

            try:
                # arcpy.env.addOutputsToMap = False
                arcpy.AddMessage("\n" + "Creating charts for output." + "\n")
                for zones in uniquezones:

                    expression = """"RaBET_ID" = '{0}'""".format(zones)
                    tview = os.path.join("in_memory", arcpy.CreateUniqueName("tv"))
                    arcpy.management.MakeTableView(outtable, tview, expression)
                    expression = """"RaBET_ID" = '{0}'""".format(zones)
                    # yearcount = str(int(arcpy.management.GetCount(tview).getOutput(0)))
                    try:
                        title = str([row[0] for row in arcpy.da.SearchCursor(tview, "RaBET_NAME", expression)][0])

                        arcpy.AddMessage("\n" + "Creating chart for: " + title)
                        # arcpy.AddMessage("yearcount")
                        # arcpy.AddMessage(yearcount)

                        graphname = outname + "_Chart_" + zones
                        graphimage = os.path.join(outdir, graphname + ".jpg")
                        maxmean = max(meanvalues)

                        if (maxmean >= 0) and (maxmean < 25):
                            graphtemplate = os.path.join(tooldatapath, "GraphTemplateMean25.tee")

                        elif (maxmean >= 25) and (maxmean < 50):
                            graphtemplate = os.path.join(tooldatapath, "GraphTemplateMean50.tee")

                        elif (maxmean >= 50) and (maxmean <= 100):
                            graphtemplate = os.path.join(tooldatapath, "GraphTemplateMean100.tee")

                        else:
                            arcpy.AddMessage("Error in Chart Data")

                        fieldID = "MEAN"

                        # Create the graph
                        graph = arcpy.Graph()

                        # Add a vertical bar series to the graph
                        graph.addSeriesBarVertical(tview, fieldID, "", "Year")

                        # Specify the title of the Graph
                        graph.graphPropsGeneral.title = title

                        # Output a graph, which is created in-memory
                        arcpy.management.MakeGraph(graphtemplate, graph, graphname)

                        # Save the graph as an image
                        if not os.path.isfile(graphimage):
                            arcpy.management.SaveGraph(graphname, graphimage, "MAINTAIN_ASPECT_RATIO", "600", "375")

                        else:
                            arcpy.AddWarning(graphimage + " already exists." + "\n" + "File not saved.")

                        # Display charts in default browser
                        if zonenum <= zonethresh:
                            webbrowser.open_new(graphimage)

                    except:
                        arcpy.AddMessage("\n" + "No Data for chart.")
                        pass
            except:
                arcpy.AddError("Error creating bar charts.")
                pass
        if zonenum > zonethresh:
            arcpy.AddMessage("\n" + "Number of polygons exceeds threshold of " + str(
                zonethresh) + "." + "\n" + "Charts will not be displayed, but were saved to output directory." + "\n" + "\n")

        arcpy.management.Delete("in_memory")

        dropFields = ["Zone_ID", "RaBET_ID", "RaBET_NUM"]
        arcpy.management.DeleteField(outlayer, dropFields)

        # Add polgons to the data frame
        try:
            mxd = arcpy.mp.MapDocument("CURRENT")
            dataFrame = mxd.listMaps("*")[0]
            addlayer = arcpy.mp.LayerFile(outlayer)
            arcpy.mp.AddLayer(dataFrame, addlayer)
            # addlayer.labelClasses[0].expression = "[RaBET_NAME]"
            # addlayer.showLabels = True
            arcpy.RefreshActiveView()
            arcpy.RefreshTOC()
        except:
            arcpy.AddWarning("Error adding layer to dataframe. Layer can be added from output directory")
            pass

        return