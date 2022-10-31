##V8 4/1/2019 - updated MLRA Data polygon and full implementation of Meta Data Layers and Layer Package output
##V7 added export of layers as layerpackage - this version only - MK
##V7 Modified by MK 10/3/2018 - Added image meta data functionality - Polygon generated with image meta data
##V6 7/20/20 Modified by MK to accept image years in UI back to 1987 - Line 99
##V6 Modified by MK 03/09/2018 added exponent term attribute table and dynamic equations.
##V5 Modified by MK 11/3/2017 added functionality  for multiple linear regression equation construction through attribute table.
##V4 Modified by MK 09/25/2017 to incorporate New Landsat SR Product as of 4/2017 -Collection 1 SR-
## Not compatible with prior Landsat formats
import arcpy
import os
import sys



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "RaBET"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [WCMapTool, RaBETAnalysisTool]


class WCMapTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Generate WC Maps"
        self.description = "This tool produces woody canopy cover maps from Landsat multispectral imagery for a given MLRA"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # First parameter
        param0 = arcpy.Parameter(
            displayName="Landsat Image Directory",
            name="imagedirin",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        # Second parameter
        param1 = arcpy.Parameter(
            displayName="MLRA Symbol",
            name="mlrain",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        # Set a value list of MLRA Names
        param1.filter.type = "ValueList"
        param1.filter.list = []
        #param1.filter.list = ['41']

        # Third parameter
        param2 = arcpy.Parameter(
            displayName="Year",
            name="yearin",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        # Set a value list of available image years
        param2.filter.type = "ValueList"
        #param2.filter.list = ['1984','1985','1986','1987','1988','1989','1990','1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2011','2012','2013','2014','2015','2016','2017']
        param2.filter.list = []

        # Fourth parameter
        param3 = arcpy.Parameter(
            displayName="Output Directory",
            name="outdir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")


        # Fifth parameter
        param4 = arcpy.Parameter(
            displayName="Area of Interest Shapefile",
            name="ROIshp",
            datatype="DEShapefile",
            parameterType="Optional",
            direction="Input")

        params = [param0, param1, param2, param3, param4]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        import os

        StartYear = 1987 # Earliest Image year selectable for input
        CompLength = 4 # Number of years used in image composition



        if parameters[0].altered:
            rootpath = sys.path[0]
            tooldatapath = os.path.join(rootpath, r"ToolData")
            mlradata = os.path.join(tooldatapath, r"MLRA_Data.shp")

            MLRAlist = [row[0] for row in arcpy.da.SearchCursor(mlradata, "MLRA_ID")]
            parameters[1].filter.list = sorted(set(MLRAlist))
            imagedir = str(parameters[0].value)
            subdirs = os.listdir(imagedir)
            subdirs = [row for row in subdirs if len(str(row)) > 30]

        if parameters[1].altered:

            year = sorted(set([x[10:14] for x in subdirs]))
            MinYear = StartYear 
            CompYears = [i for i in year if int(i) >= MinYear]


            parameters[2].filter.list = CompYears

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""



        # Import modules
        import os
        import time
        import sys
        import arcpy
        import glob
        import datetime
        import shutil
        from operator import itemgetter
        from arcpy import env
        from arcpy.sa import *

        # Version and update information
        version = "RaBET 6.0 demo"
        modified = r"03/09/2018"
        executiontime = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")

        # Define input parameters
        imagedir = parameters[0].valueAsText # Landsat image directory
        mlraID = parameters[1].valueAsText # MLRA name/ID
        yearID = str(parameters[2].valueAsText) # Year to process
        outdir = parameters[3].valueAsText # Qutput directory
        ROIshp = parameters[4].valueAsText # Option ROI shapefile

        CompLength = 4 # Number of years in the vegetation index composite

        startyear = str(int(yearID) - (CompLength -1))
        endyear = yearID






        # Define additional parameters
        bufferdist = "10 Kilometers" # Buffer around ROI for cloud cover check (include units)

        # Check out the Spatial Analyst extension
        try:
            arcpy.CheckOutExtension("Spatial")
        except:
            arcpy.AddWarning("Spatial Analysis licence is not available.")
            pass

        #Determine directories containg Python toolbox and reource data
        rootpath = sys.path[0]
        tooldatapath = os.path.join(rootpath, r"ToolData")
        mlradata = os.path.join(tooldatapath, r"MLRA_Data.shp")
        scratchdirname = "scratch"
        scratchwsname = "scratch.gdb"
        scratchdir = os.path.join(rootpath, scratchdirname)
        scratchws = os.path.join(rootpath, scratchdirname, scratchwsname)
        
        NLCDmaskloc = os.path.join(tooldatapath, "NLCD_Mask_" + mlraID + ".tif")
        StateMask84B = os.path.join(tooldatapath, "State_Mask_84B.tif")
        StateMask85 = os.path.join(tooldatapath, "State_Mask_85.tif")
        WCarchivedirname = "WC_Image_Archive"
        WCarchivegdbname = "WC_Image_Archive.gdb"
        WCarchivedir = os.path.join(rootpath,WCarchivedirname)
        WCarchivegdb = os.path.join(rootpath,WCarchivedirname,WCarchivegdbname)
        metadir = os.path.join(rootpath, WCarchivedirname, "metadata")
        WCimagename = "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear
        WCarchiveimage = os.path.join(WCarchivegdb, WCimagename)


        # Set Geoprocessing environments
        arcpy.env.scratchWorkspace = scratchws
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = outdir

        #Display general processing information to the console
        arcpy.AddMessage("\n" + "NRCS ArcGIS Woody Cover Tool Version " + version + "\n" + "Last modified on " + modified + "\n")
        arcpy.AddMessage("Output Directory: " + "\n" + outdir + "\n")
        arcpy.AddMessage("Processing images for:")
        arcpy.AddMessage("MLRA: " + mlraID)
        arcpy.AddMessage("Year: " + yearID +  "\n")
        arcpy.AddMessage("Composite using years: " + startyear + " - " + endyear + "." + "\n")
        arcpy.AddMessage(WCimagename)


        # Check if scratch dir and gdb exist. Clear existing data in gdb.
        if not arcpy.Exists(scratchws):
            arcpy.AddMessage("Scratch GDB DNE")
            if not arcpy.Exists(scratchdir):
                arcpy.AddMessage("Scratch Dir DNE")
                arcpy.CreateFolder_management(rootpath, scratchdirname)
            arcpy.CreateFileGDB_management(scratchdirname, scratchwsname, "10.0")

        else:
            arcpy.AddMessage("Scratch GDB exists: "+ scratchws)
            scratchfiles = os.listdir(scratchws)
            
            try:
                
                arcpy.Delete_management(scratchws)
                arcpy.CreateFileGDB_management(scratchdirname, scratchwsname, "10.0")
                
                arcpy.AddMessage("Scratch directory recreated.")
            except:
                arcpy.AddMessage("Scratch directory could not be deleted.")




        # Check if image archive exists. Check if image has been processed.
        if not arcpy.Exists(WCarchivegdb):
            arcpy.AddMessage("WC GDB DNE")
            if not arcpy.Exists(WCarchivedir):
                arcpy.AddMessage("WC Dir DNE")
                arcpy.CreateFolder_management(rootpath, WCarchivedirname)
            arcpy.CreateFileGDB_management(WCarchivedir, WCarchivegdbname, "10.0")

        else:
            arcpy.AddMessage("Image archive exists")

        if not arcpy.Exists(metadir):
            arcpy.CreateFolder_management(WCarchivedir, "metadata")

        if not arcpy.Exists(WCarchiveimage):
            arcpy.AddMessage(WCimagename + " does not exist in archive. Image will be processed.")







          ####################### MSAVI PROCESSING ##########################



            # Check for region of interest and determine tiles to be processed

            if ROIshp is None:
                ROIexist = False
                arcpy.AddMessage("No region of interest shapefile was defined." + "\n")

                theregion = mlradata
                regionname = "MLRA " + mlraID

                metaname = os.path.join(outdir, "meta_" + mlraID + "_" + str(yearID) + ".txt")
            else:
                ROIexist = True
                ROIname = os.path.splitext(os.path.basename(ROIshp))[0]
                arcpy.AddMessage("Subsetting region of interest from shapefile: " + "\n" + str(ROIshp) + "\n")

                inttest = os.path.join("in_memory", "inttest")
                arcpy.Intersect_analysis([mlradata, ROIshp], inttest)
                rowcount = int(arcpy.GetCount_management(inttest).getOutput(0))

                metaname = os.path.join(outdir, "meta_" + mlraID + "_" + str(yearID) + "_" + ROIname + ".txt")

                if rowcount == 0:
                    arcpy.AddError("Error subsetting MLRA " + mlraID + "." + "\n" + "Ensure that the area of interest shapefile intersects the MLRA.")
                    sys.exit(0)
                else:
                    # Clip MLRA_Data.shp to region of interest
                    mlradataROI = os.path.join("in_memory", "MLRA_Data_ROI")
                    arcpy.Clip_analysis(mlradata, ROIshp, mlradataROI)
                    theregion = mlradataROI
                    regionname = ROIname

            # List tiles containing region of interest
            expression = """"MLRA_ID" = '{0}'""".format(mlraID)
            pathrowlist = [row[0] for row in arcpy.da.SearchCursor(theregion, "Path_Row", expression)]

            # List start and end DOY constraints
            DOYstartlist = [row[0] for row in arcpy.da.SearchCursor(theregion, "DOY_Start", expression)]
            DOYendlist = [row[0] for row in arcpy.da.SearchCursor(theregion, "DOY_End", expression)]

            # List Landsat subdirectories containing surface imagery
            subdirs = os.listdir(imagedir)
            subdirs = [row for row in subdirs if len(str(row)) > 30] # Check if filename is at least 30 characters long

            # Parse variables from directory strings
            pathrows = [x[4:10] for x in subdirs]
            years = [x[10:14] for x in subdirs]
            
            months = [x[14:16] for x in subdirs]
            days = [x[16:18] for x in subdirs]



            #Calculate DOY from day month year
            DOYs = []
            for yyyy, mm, dd in zip(years, months, days):
                yyyy = int(yyyy)
                mm = int(mm)
                dd = int(dd)
                DayofYear = datetime.date(yyyy, mm, dd).timetuple().tm_yday
                DOYs.append(DayofYear)

            

            

            

            # Screen image directories for path/row
            pathrowindex = []
            for tile in pathrowlist:
                pathrowindex = pathrowindex + [row for row, dirs in enumerate(pathrows) if dirs == str(tile)]

            # Screen image directoies for year
            yearindex = [row for row, dirs in enumerate(years) if (int(dirs) >= int(startyear)) and (int(dirs) <= int(endyear))]

            # Screen image directories for day of year range
            DOYindex = []
            for startday,endday in zip(DOYstartlist, DOYendlist):
                DOYindex = DOYindex + [row for row, dirs in enumerate(DOYs) if (int(dirs) >= int(startday)) and (int(dirs) <= int(endday))]

            # Determine available images based on intersection of image constraint indices (year, P/R, DOY)
            indexintersection = set.intersection(set(pathrowindex), set(yearindex), set(DOYindex))
            availableimages = [subdirs[row] for row in indexintersection]
            availableimages.sort()
            availableimagenum = len(availableimages)

            # Set of Unique Landsat tiles with available images in region of interest
            uniquetiles = sorted(set([subdirs[row][4:10] for row in indexintersection]))
            numbertiles = len(uniquetiles)

#create metadata text file
            
            metadata = open(str(metaname),"w")
            metadata.write("Version: " + version + "\n")
            metadata.write("Image created at " + executiontime + "\n"  )
            metadata.write("\n")
            metadata.write("RaBET user inputs: " + "\n")
            metadata.write("    Landsat image directory: " + imagedir  + "\n")
            metadata.write("    MLRA symbol: " + mlraID  + "\n")
            metadata.write("    Input year: " + yearID  + "\n")
            metadata.write("    Output directory: " + outdir  + "\n")
            if ROIexist == True:
                metadata.write("    Subset to region of interest: " + ROIname  + "\n")

            metadata.write("\n")
            metadata.write("Image composite settings:" + "\n")
            metadata.write("    Composite length: " + str(CompLength)  + "\n")
            metadata.write("    Composite start year: " + str(startyear)  + "\n")
            metadata.write("    Composite end year: " + str(endyear)  + "\n")
            metadata.write("\n")


            # Display image constraint information to the console.
            if availableimagenum > 0:
                arcpy.AddMessage("Landsat WRS-2 Path/Row IDs in " + regionname  +  ": " + ", ".join(map(str, uniquetiles)))
                arcpy.AddMessage("The following " + str(availableimagenum) + " image(s) match date and location criteria:")
                arcpy.AddMessage("\n".join(availableimages) + "\n")
                metadata.write("Landsat WRS-2 Path/Row IDs in " + regionname  +  ": " + ", ".join(map(str, uniquetiles)) + "\n")
                metadata.write("\n")
                metadata.write("The following " + str(availableimagenum) + " image(s) match date and location criteria:" + "\n")
                metadata.write("\n".join(availableimages) + "\n")
                metadata.write("\n")
                metadata.write("*******************************************************************************")
                metadata.write("\n")
                
            else:
                arcpy.AddWarning("No valid Landsat images avalailable for " + yearID + "." + "\n")
                metadata.write("No valid Landsat images avalailable for " + yearID + "." + "\n")
                sys.exit(0)


            # Landsat image geoprocessing
            pathrowprocessed = [] #List of proceesed tiles
            finalimages = []
            inmemorylist = []


            

        #imagedir = parameters[0].valueAsText # Landsat image directory
        #mlraID = parameters[1].valueAsText # MLRA name/ID
        #yearID = str(parameters[2].valueAsText) # Year to process
        #outdir = parameters[3].valueAsText # Qutput directory
        #ROIshp = parameters[4].valueAsText # Option ROI shapefile

        #CompLength = 4 # Number of years in the vegetation index composite

        #startyear = str(int(yearID) - (CompLength -1))
        #endyear = yearID
            
            


            for tiles in uniquetiles:
                # Create list to store dictionaries of image properties
                filteredimages = []
                tileimageslist = []
                metadata.write("\n" + "Composite images for path/row ID: " + tiles + ":" + "\n")


                arcpy.AddMessage("\n" + r"Processing path/row ID: " + tiles + "\n")

                for imindex, images in enumerate(availableimages):

                    # start time for image processing
                    totalstarttime = time.time()

                    # Define variables from Landsat directory string
                    currentpathrow = images[4:10] # path/row of current image in loop
                    currentimagename = images[0:20]
                    satellite = images[3:4]
                    pathname = images[5:7]
                    rowname = images[8:10]
                    month = images[14:16]
                    day = images[16:18]
                    year = images[10:14]
                   
                    DOY = str(datetime.date(int(year), int(month), int(day)).timetuple().tm_yday)


                    

                    while currentpathrow == tiles:

                        # Define current Landsat image directory
                        imagenum = imindex + 1
                        arcpy.AddMessage("\n" + "Processing image " + str(imagenum) + " of " + str(availableimagenum) + ": " + images + "\n")
                        currentimagedir = os.path.join(imagedir, images)

                        # Cloud Cover Screening

                        processstarttime = time.time()

                        # Load cloud mask image from current directory as raster
                        cfmask = glob.glob(currentimagedir  + '\\*pixel_qa.tif')
                        cloudraster =  arcpy.Raster(cfmask[0])

                        # Define projection to current cloud mask image
                        arcpy.env.outputCoordinateSystem = cloudraster
						# Updated 2021/04/02
                        # Set the snap raster for the processing because otherwise the Clip tool will shift the pixels
                        # when the clipping extent from the extentcloud and clipcloud inputs is different than 
                        # the extent of the dataset to be clipped
                        arcpy.env.snapRaster = cloudraster


                        # Define Pixel QA Band Mask Definitions - Includes Landsat 5 & 8 Defs
                        
                        clearmap = [[66,0], [130,0], [322,0], [386,0]] # Class 0 - clear
                        
                        watermap = [[68,1], [132,1], [324,1], [388,1]] # Class 1 - water
                        
                        cloudshadowmap = [[72,2], [136,2], [328,2], [392,2]] # Class 2 - shadow
                        
                        snowmap = [[80,3], [144,3], [336,3], [400,3]] # Class 3 - snow
                        
                        cloudmap = [[96,4], [112,4], [160,4], [176,4], [224,4], [352,4], [368,4], [416,4], [432,4], [480,4], [832,4], [836,4], [840,4], [848,4], [864,4], [880,4], [900,4], [904,4], [912,4], [928,4], [944,4], [992,4]] # Class 4 - cloud


                        # Reclass map of all class changes
                        classmap = clearmap + watermap + cloudshadowmap + snowmap + cloudmap
                        remap = RemapValue(classmap)

                        #Reclass Pixel QA image to 0-4 classes
                        cloudrasterRecl = Reclassify(cloudraster, "Value", remap, "DATA")
                        QAreclLoc = os.path.join(scratchws, currentimagename + "_cfmask_recl")
                        cloudrasterRecl.save(QAreclLoc)

                        pathrowprocessed.append(currentpathrow)

                        # Select polygon for current path/row
                        expression = """"Path_Row" = '{0}' AND "MLRA_ID" = '{1}'""".format(currentpathrow, mlraID)
                        cliplayertile = "Clip_Layer_" + currentpathrow
                        arcpy.MakeFeatureLayer_management(mlradata, cliplayertile, expression)

                        #Define geometry/extent parameters for region of analysis
                        clipenvelopetile = os.path.join("in_memory", "Clip_Envelope_" + currentpathrow )
                        arcpy.FeatureEnvelopeToPolygon_management(cliplayertile, clipenvelopetile, "SINGLEPART")
                        desccliptile = arcpy.Describe(clipenvelopetile)
                        extentcliptile = str(desccliptile.extent).translate(None, 'NaN')

                        #Define geometry parameters for subset region if it exists
                        if ROIexist:

                            cliplayerROI = arcpy.CreateUniqueName(os.path.join("in_memory", "Clip_Layer_" + currentpathrow))
                            arcpy.MakeFeatureLayer_management(mlradataROI, cliplayerROI, expression)

                            clipenvelopeROI = os.path.join("in_memory", "Clip_Envelope_" + currentpathrow )
                            arcpy.FeatureEnvelopeToPolygon_management(cliplayerROI, clipenvelopeROI, "SINGLEPART")
                            descclipROI = arcpy.Describe(clipenvelopeROI)
                            extentclipROI = str(descclipROI.extent).translate(None, 'NaN')

                            cloudbuffer = os.path.join("in_memory", "Cloud_Buffer_" + currentpathrow )
                            arcpy.Buffer_analysis(cliplayerROI, cloudbuffer, bufferdist, "#", "#", "All", "MLRA_ID")

                            cloudbufferint = os.path.join("in_memory", "Cloud_Buffer_Intersect_" + currentpathrow )
                            arcpy.Clip_analysis(cloudbuffer, cliplayertile, cloudbufferint)

                            bufferenvelope = os.path.join("in_memory", "Buffer_Envelope_" + currentpathrow)
                            arcpy.FeatureEnvelopeToPolygon_management(cloudbufferint, bufferenvelope, "SINGLEPART")

                            descbuffer = arcpy.Describe(bufferenvelope)
                            extentbuffer = str(descbuffer.extent).translate(None, 'NaN')

                            extentcloud = extentbuffer
                            clipcloud = cloudbufferint

                            extentVI = extentclipROI
                            clipVI = cliplayerROI

                        else:
                            extentcloud =  extentcliptile
                            clipcloud = cliplayertile

                            extentVI = extentcliptile
                            clipVI = cliplayertile

                        # Begin cloud screening using cfmask Landsat rasters
                        arcpy.AddMessage(r"Screening image for cloud/snow cover.")
                        cloudrasterclip = os.path.join(scratchws, currentimagename + "_cfmask_clip")

                        arcpy.Clip_management(cloudrasterRecl, extentcloud, cloudrasterclip, clipcloud,  "255", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
                        arcpy.BuildRasterAttributeTable_management(cloudrasterclip, "Overwrite")


                        cursorcloud = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count")]
                        cursorvalues = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Value")]

                        totalpixels = sum(cursorcloud)



                        # itemize CFmask classes within region of interest
                        
                        if 0 in cursorvalues:
                            expression = """"Value" = 0"""
                            clearpixels = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count", expression)]

                        else:
                            clearpixels = [0]

                        if 1 in cursorvalues:
                            expression = """"Value" = 1"""
                            waterpixels = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count", expression)]

                        else:
                            waterpixels = [0]

                        if 2 in cursorvalues:
                            expression = """"Value" = 2"""
                            cloudshadowpixels = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count", expression)]

                        else:
                            cloudshadowpixels = [0]

                        if 3 in cursorvalues:
                            expression = """"Value" = 3"""
                            snowpixels = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count", expression)]

                        else:
                            snowpixels = [0]

                        if 4 in cursorvalues:
                            expression = """"Value" = 4"""
                            cloudpixels = [row[0] for row in arcpy.da.SearchCursor(cloudrasterclip, "Count", expression)]

                        else:
                            cloudpixels = [0]

                        # Calculate CFmask percent coverage for all classes
                        clearpercent = 100 * clearpixels[0] / totalpixels
                        waterpercent = 100 * waterpixels[0] / totalpixels
                        cloudshadowpercent = 100 * cloudshadowpixels[0] / totalpixels
                        snowpercent = 100 * snowpixels[0] / totalpixels
                        cloudpercent = 100 * cloudpixels[0] / totalpixels
                        totaldisturbedpercent = cloudpercent + cloudshadowpercent + snowpercent

                    
                        arcpy.AddMessage(str(int(round(cloudpercent))) + "% cloud pixels")
                        arcpy.AddMessage(str(int(round(cloudshadowpercent))) + "% cloud shadow pixels")
                        arcpy.AddMessage(str(int(round(snowpercent))) + "% snow pixels" + "\n")

                        # Determine cloud cover threshold for tile
                        expression = """"Path_Row" = '{0}' AND "MLRA_ID" = '{1}'""".format(currentpathrow, mlraID)
                        cloudthresh = [row[0] for row in arcpy.da.SearchCursor(theregion, "Cloud_Thr", expression)]
                        
                        

                        



                        processstarttime = time.time()
                        processendtime = time.time() - processstarttime
                        processendtime = int(round(processendtime))
                        arcpy.Delete_management(cloudrasterclip)

                        

                        # check to see if image meets cloud cover threshold criteria
                        if int(totaldisturbedpercent) > int(cloudthresh[0]):
                            arcpy.AddMessage("Total disturbed pixels (snow, cloud and shadow) = " + str(int(round(totaldisturbedpercent))) + "%. " + "\n" + "Image exceeds cloud cover threshold of " + str(cloudthresh[0]) + r"%." )
                            metadata.write("    " + images + " - Excluded - Total disturbed pixels: " + str(int(round(totaldisturbedpercent))) + "%. " +  "\n")

                            

                        # Begin vegetation index processing
                        else:
                            arcpy.AddMessage("Total disturbed pixels (snow, cloud and shadow) = " + str(int(round(totaldisturbedpercent))) + "%. " + "\n" + "Image meets cloud cover criteria (" + str(cloudthresh[0]) + r"%)." )
                            processstarttime = time.time()
                            
                            
                            #Get number of subregions in MLRA
                            SubregionsNum = int([row[0] for row in arcpy.da.SearchCursor(theregion, "Subregions", expression)][0])

                            
                            indexlist = []
                            for RegionNum in range(1,SubregionsNum + 1):
                                indexstring = [row[0] for row in arcpy.da.SearchCursor(theregion, "Index_" + str(RegionNum), expression)][0]
                                indexlist.extend([x.strip() for x in indexstring.split(',')])

                            indexlist = set(indexlist)
                            indexnumber = len(indexlist)


                            arcpy.AddMessage("\n" + "Processing " + str(indexnumber) + " vegetation indices: " + indexstring + ". " )
                            for indexname in indexlist:

                                # select appropriate Landsat bands and apply conversion factor
                                if satellite == '8':
                                    blueimage = glob.glob(currentimagedir + '\\*sr_band2.tif')
                                    greenimage = glob.glob(currentimagedir + '\\*sr_band3.tif')
                                    redimage = glob.glob(currentimagedir + '\\*sr_band4.tif')
                                    nirimage = glob.glob(currentimagedir + '\\*sr_band5.tif')
                                    swir1image = glob.glob(currentimagedir + '\\*sr_band6.tif')
                                    swir2image = glob.glob(currentimagedir + '\\*sr_band7.tif')
                                    scalefactor = 0.0001

                                elif satellite == '5' or satellite == '7':
                                    blueimage = glob.glob(currentimagedir + '\\*sr_band1.tif')
                                    greenimage = glob.glob(currentimagedir + '\\*sr_band2.tif')
                                    redimage = glob.glob(currentimagedir + '\\*sr_band3.tif')
                                    nirimage = glob.glob(currentimagedir + '\\*sr_band4.tif')
                                    swir1image = glob.glob(currentimagedir + '\\*sr_band5.tif')
                                    swir2image = glob.glob(currentimagedir + '\\*sr_band7.tif')
                                    scalefactor = 0.0001

                                else:
                                    arcpy.AddMessage('Unknown Satellite Platform')

                                
                                arcpy.AddMessage("Calculating " + indexname + " for Landsat " + str(satellite) + " image " + currentimagename + "." )


                                # Calculate vegetation index
                                if indexname == "MSAVI":
                                    redband = arcpy.Raster(redimage[0])
                                    nirband = arcpy.Raster(nirimage[0])
                                    NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    VI = Con(cloudrasterRecl == 0, ((((2 * (scalefactor * nirband)) + 1) -  (SquareRoot( Square(2 * (scalefactor * nirband) + 1) - 8 * ((scalefactor * nirband) - (scalefactor * redband))))) / 2))
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    VImask = Con(NLCDmask == 1, VI)
                                elif indexname == "GSATVI":
                                    greenband = arcpy.Raster(greenimage[0])
                                    nirband = arcpy.Raster(nirimage[0])
                                    swir1band = arcpy.Raster(swir1image[0])
                                    #Get soil adjustment factor L from MLRA Data
                                    L = float([row[0] for row in arcpy.da.SearchCursor(theregion, "L", expression)][0])
                                    arcpy.AddMessage("Soil adjustment factor, L: " + str(L) + "." )
                                    # L = 0.1
                                    NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    VI = Con(cloudrasterRecl == 0, (((((scalefactor * nirband) - (scalefactor * greenband)) / ((scalefactor * nirband) + (scalefactor * greenband) + L)) * (1 + L)) - ((scalefactor * swir1band) / 2)))
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    VImask = Con(NLCDmask == 1, VI)
                                elif indexname == "NDI5":
                                    nirband = arcpy.Raster(nirimage[0])
                                    swir1band = arcpy.Raster(swir1image[0])
                                    NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    VI = Con(cloudrasterRecl == 0, ((scalefactor * swir1band) - (scalefactor * nirband)) / ((scalefactor * swir1band) + (scalefactor * nirband)))
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    VImask = Con(NLCDmask == 1, VI)
                                elif indexname == "NDVI":
                                    redband = arcpy.Raster(redimage[0])
                                    nirband = arcpy.Raster(nirimage[0])
                                    NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    VI = Con(cloudrasterRecl == 0, ((scalefactor * nirband) - (scalefactor * redband)) / ((scalefactor * nirband) + (scalefactor * redband)))
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    VImask = Con(NLCDmask == 1, VI)
                                elif indexname == "Constant":
                                    #redband = arcpy.Raster(redimage[0])
                                    #NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    #VI = Con(cloudraster == 0, redband / redband)
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    #VImask = Con(NLCDmask == 1, VI)
                                    redband = arcpy.Raster(redimage[0])
                                    
                                    NLCDmask = arcpy.Raster(NLCDmaskloc)
                                    VI = Con(cloudrasterRecl == 0, (redband) / (redband) )
                                    # Apply urban, barren, water, ag mask from NLCD image
                                    VImask = Con(NLCDmask == 1, VI)
                                
                                else:
                                    arcpy.AddMessage("Vegetation index does not exist")

                                # Calculate zonal statistics
                                zonalstats = os.path.join("in_memory", "z" + pathname + rowname + DOY )
                                ZonalStatisticsAsTable(cliplayertile, "ZoneID", VImask, zonalstats, "DATA", "ALL")
                                meanvi = str([row[0] for row in arcpy.da.SearchCursor(zonalstats, "MEAN")][0])
                                stdvi = [row[0] for row in arcpy.da.SearchCursor(zonalstats, "STD")]
                                area = [row[0] for row in arcpy.da.SearchCursor(zonalstats, "AREA")]
                                arcpy.Delete_management(zonalstats)

                                #Clip VI raster to region of interest
                                VIcliptile = os.path.join(scratchws, indexname + pathname + rowname + DOY)
                                arcpy.Clip_management(VImask, extentVI, VIcliptile, clipVI, "255", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
                                inmemorylist.append(VIcliptile)



                                wcname = os.path.join(scratchws, "wc" + mlraID + pathname + rowname)
                                


                                # Add image properties to dictionary
                                imageproperties = {
                                    "current image": currentimagename,
                                    "path": pathname,
                                    "row": rowname,
                                    "DOY": DOY,
                                    "year": year,
                                    "vegetation index": indexname,
                                    "mean VI": meanvi,
                                    "stddv VI": stdvi,
                                    "area": area,
                                    "Cliped VI": VIcliptile,
                                    "Cloud Percent": cloudpercent,
                                    "MLRA": mlraID,
                                    "WC Image": wcname
                                    }
                                filteredimages.append(imageproperties)
                                tileimageslist.append(VIcliptile)
                                

                                pathnamestr = pathname
                                rownamestr = rowname
                                pathrowname = "0" + pathname + "0" + rowname


                                processendtime = time.time() - processstarttime
                                processendtime = int(round(processendtime))
                                

                            metadata.write("    " + images + " - Total disturbed pixels (snow, cloud and shadow) = " + str(int(round(totaldisturbedpercent))) + "%. " +  "\n")
                            arcpy.Delete_management(cloudrasterRecl)
                                
                        break

                # Calculate woody canopy cover if images exist
                if filteredimages:
                    processstarttime = time.time()

                    cellstattype = [row[0] for row in arcpy.da.SearchCursor(theregion, "Cell_Stat", expression)][0]
                    cellstatsmetric = str(cellstattype)
                    
                    compositelist = []
                    arcpy.AddMessage("\n" + "Compositing images for " + pathnamestr + "/" + rownamestr + ".")


                    for indexname in indexlist:
                        arcpy.AddMessage("Compositing " + indexname + " images"  + " using " + cellstatsmetric + "." )
                        metadata.write("\n" + indexname + " composite method: "  + "\n")
                        metadata.write("    " +  cellstatsmetric + "\n")
                        
                        indeximageslist = []               
                        for prtiles in tileimageslist:
                            if indexname in prtiles:
                                indeximageslist.append(prtiles)
                                
                      
                        
                        compositename = os.path.join(scratchws, indexname + pathnamestr + rownamestr + "CMP")
                        outCellStats = CellStatistics(indeximageslist, cellstatsmetric, "DATA")
                        outCellStats.save(compositename)
                        compositelist.append(compositename)
                        inmemorylist.append(compositename)





                        


                    arcpy.AddMessage("\n" + "Calculating woody canopy cover for MLRA " + mlraID + "-" + pathnamestr + "/" + rownamestr + "." + "\n")

                    if SubregionsNum > 3:
                        arcpy.AddError("Number of subregions exceeds maximum of 3.")
                        sys.exit(0)
                    
                    if SubregionsNum > 1:
                        try:
                            RasterClass = Raster(os.path.join(tooldatapath, "Subregions"  + "_" + mlraID + ".tif"))

                        except:
                            arcpy.AddWarning("Subregions raster does not exist.")
                            pass

                    for CalcRegionNum in range(1,SubregionsNum + 1):
                        CalcCoeffString = [row[0] for row in arcpy.da.SearchCursor(theregion, "Coeff_" + str(CalcRegionNum), expression)][0]
                        CalcCoeffList = [float(x.replace(" ", "")) for x in CalcCoeffString.split(',')]
                        CalcCoeffStrList = [x.replace(" ", "") for x in CalcCoeffString.split(',')]
                        CalcCoeffNum = len(CalcCoeffList)
                        CalcIndexString = [row[0] for row in arcpy.da.SearchCursor(theregion, "Index_" + str(CalcRegionNum), expression)][0]
                        CalcIndexList = [x.strip() for x in CalcIndexString.split(',')]
                        CalcIndexNum = len(CalcIndexList)
                        CalcExpString = [row[0] for row in arcpy.da.SearchCursor(theregion, "Exp_" + str(CalcRegionNum), expression)][0]
                        CalcExpList = [float(x.replace(" ", "")) for x in CalcExpString.split(',')]
                        CalcExpStrList = [x.replace(" ", "") for x in CalcExpString.split(',')]
                        CalcExpNum = len(CalcExpList)


                        if CalcCoeffNum > CalcIndexNum:
                            arcpy.AddError("Number of vegetation indices does not match the number of coefficients.")
                            sys.exit(0)

                        if CalcCoeffNum > 4:
                            arcpy.AddError("Number of tems in WC equation exceeds maximum of 3.")
                            sys.exit(0)


                        CalcCompositeList = []
                        
                        
                        for IndexCount, VegIndex in enumerate(CalcIndexList):
                            CalcCompositeList.extend([comp for comp in compositelist if CalcIndexList[IndexCount] in comp])

                            

                        if CalcRegionNum == 1:
                            if CalcCoeffNum == 1:
                                arcpy.AddMessage("Subregion 1")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])

                                Region1 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")

                            if CalcCoeffNum == 2:
                                arcpy.AddMessage("Subregion 1")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])

                                Region1 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")

                            if CalcCoeffNum == 3:
                                arcpy.AddMessage("Subregion 1")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])

                                Region1 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                            if CalcCoeffNum == 4:
                                arcpy.AddMessage("Subregion 1")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + "  + CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])
                                Raster4 = Raster(CalcCompositeList[3])

                                Region1 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3  ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + " +  CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")
                                metadata.write("    GSATVI L-factor = " + str(L) + "\n")

                        if CalcRegionNum == 2:
                            if CalcCoeffNum == 1:
                                arcpy.AddMessage("Subregion 2")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])

                                Region2 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")


                            if CalcCoeffNum == 2:
                                arcpy.AddMessage("Subregion 2")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])

                                Region2 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")


                            if CalcCoeffNum == 3:
                                arcpy.AddMessage("Subregion 2")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])

                                Region2 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                            if CalcCoeffNum == 4:
                                arcpy.AddMessage("Subregion 2")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + "  + CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])
                                Raster4 = Raster(CalcCompositeList[3])

                                Region2 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3  ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + " +  CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")
                                metadata.write("    GSATVI L-factor = " + str(L) + "\n")

                        if CalcRegionNum == 3:
                            if CalcCoeffNum == 1:
                                arcpy.AddMessage("Subregion 3")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])

                                Region3 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + "\n")

                            if CalcCoeffNum == 2:
                                arcpy.AddMessage("Subregion 3")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])

                                Region3 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]))))

                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  "\n")


                            if CalcCoeffNum == 3:
                                arcpy.AddMessage("Subregion 3")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])

                                Region3 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] +  " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] +  " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] +   "\n")

                            if CalcCoeffNum == 4:
                                arcpy.AddMessage("Subregion 3")
                                arcpy.AddMessage("WC = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + "  + CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")

                                Raster1 = Raster(CalcCompositeList[0])
                                Raster2 = Raster(CalcCompositeList[1])
                                Raster3 = Raster(CalcCompositeList[2])
                                Raster4 = Raster(CalcCompositeList[3])

                                Region3 = Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))< 0, 0, Con(Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3 ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3])) > 100, 100, Int((CalcCoeffList[0] * Raster1 ** CalcExpList[0]) + (CalcCoeffList[1] * Raster2 ** CalcExpList[1]) + (CalcCoeffList[2] * Raster3  ** CalcExpList[2]) + (CalcCoeffList[3] * Raster4 ** CalcExpList[3]))))
                            
                                metadata.write("\n" + "WC Eqn. for MLRA " + mlraID + " subregion " + str(CalcRegionNum) + ":" +  "\n")
                                metadata.write("    WC% = " + CalcCoeffStrList[0] + "*" + CalcIndexList[0] + "^" + CalcExpStrList[0] + " + "  + CalcCoeffStrList[1] + "*" + CalcIndexList[1] + "^" + CalcExpStrList[1] + " + "  + CalcCoeffStrList[2] + "*" + CalcIndexList[2] + "^" + CalcExpStrList[2] + " + " +  CalcCoeffStrList[3] + "*" + CalcIndexList[3] + "^" + CalcExpStrList[3] + "\n")
                                metadata.write("    GSATVI L-factor = " + str(L) + "\n")


                    if SubregionsNum == 1:
                        wcbounds = Region1
                    elif SubregionsNum == 2:
                        wcbounds = Con(RasterClass == 1,Region1, Con(RasterClass == 2,Region2))
                    elif SubregionsNum == 3:
                        wcbounds = Con(RasterClass == 1,Region1, Con(RasterClass == 2,Region2, Con(RasterClass == 3, Region3)))
                    
                    metadata.write("\n")
                    metadata.write("*******************************************************************************")
                    metadata.write("\n")
                    
                    # Updated 2021/04/02
                    #wcbounds.save(wcname)
                    arcpy.CopyRaster_management(in_raster=wcbounds, out_rasterdataset=wcname, pixel_type="8_BIT_UNSIGNED")
                    inmemorylist.append(wcname)

                    

                    mosaicdict = {"area": area, "WC Image": wcname, "Tile":tiles}
                                    
                    finalimages.append(mosaicdict)
                    

                    

                    processendtime = time.time() - processstarttime
                    processendtime = int(round(processendtime))

                    arcpy.AddMessage("Woody canopy cover calculation complete for MLRA " + mlraID + r", parth/row " + pathnamestr + r"/" + rownamestr + " in " + str(processendtime) + " seconds)." )

                    
                else:
                        arcpy.AddMessage("No image for the curent path/row exists.")

            arcpy.AddMessage("\n" + "Creating mosaic image for MLRA " + mlraID + "." + "\n")
            if len(finalimages) < 1:
                            arcpy.AddError("No images available for mosaic.")
                            sys.exit(0)


            #Sort by largest area
            sortedimagesarea = sorted(finalimages, key=itemgetter("area"), reverse=True)
            processedtiles = []
            
            
            finalimagenames = []
            finaltilenames =[]
            
            for dict in sortedimagesarea:
                selectedimg = dict.get("WC Image")
                selectedtile = dict.get("Tile")
                finalimagenames.append(selectedimg)
                finaltilenames.append(selectedtile)

            # Updated 2021/04/02
            # Set the output coordinate system of the final mosaic to the projection of the image that has the largest area in the MLRA
            arcpy.env.outputCoordinateSystem = arcpy.Describe(finalimagenames[0]).spatialReference


            if len(pathrowprocessed) < int(numbertiles):
                processedtilesset = set(processedtiles)
                notprocessed = [tile for tile in uniquetiles if tile not in processedtilesset]

                arcpy.AddMessage("No valid images found for the following tiles:")
                for images in notprocessed:
                    arcpy.AddMessage(images)

            #Mosaic selected images

            if ROIexist:
                wcmosaic = "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "_" + ROIname + ".tif"
                wclayertemp = os.path.join("WC_MLRA_" + mlraID + "_" + startyear + "-" + endyear + "_" + ROIname + ".lyr")
                wclayer = os.path.join(outdir, "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "_" + ROIname + ".lyr")
            else:
                wcmosaic =  "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + ".tif"
                wclayertemp = os.path.join(outdir, "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "tmp.lyr")
                wclayer = os.path.join(outdir, "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + ".lyr")
            try:
				# Updated 2021/04/02
                #arcpy.CreateRasterDataset_management(outdir, wcmosaic, "", "8_BIT_UNSIGNED", "","1", "","","","")
                inputstring = ';'.join(row for row in finalimagenames)
                tilestring = ', '.join(row for row in finaltilenames)
                targetstring =  os.path.join(outdir, wcmosaic)
            except:
                arcpy.AddWarning("Error creating mosaic. Ensure image does not exist in output directory.")

                pass

            mosaicmethod = "FIRST"

			# Updated 2021/04/02
            #arcpy.Mosaic_management(inputstring,targetstring,mosaicmethod,mosaicmethod,"", "255", "", "", "")
            arcpy.MosaicToNewRaster_management(input_rasters=inputstring, output_location=outdir, raster_dataset_name_with_extension=wcmosaic,pixel_type="8_BIT_UNSIGNED", mosaic_method=mosaicmethod, number_of_bands=1)
            arcpy.AddMessage("Mosaic created using method: " + mosaicmethod)
            arcpy.AddMessage("Tile input order: " + tilestring)

            

            metadata.write("Mosaic created using method: " + mosaicmethod + "\n")
            metadata.write("Tile input order: " + tilestring + "\n")
            if ROIexist == True:
                metadata.write("\n")
                metadata.write("*******************************************************************************" )
                metadata.write("\n")
                metadata.write("Image was subset from an archived image using: " + ROIname)
            metadata.seek(0,0)
            metadata.close

            if ROIshp is None:
                arcpy.RasterToGeodatabase_conversion(targetstring, WCarchivegdb)
                metanamearchive = os.path.join(metadir, "meta_" + mlraID + "_" + str(yearID) + ".txt")
                metaname = os.path.join(outdir, "meta_" + mlraID + "_" + str(yearID) + ".txt")
                shutil.copy(metaname, metanamearchive)

            #Add layer and symbology
            arcpy.MakeRasterLayer_management(targetstring, wclayertemp, "#", "#", "1")
            wcsymbology = os.path.join(tooldatapath, "wcsymbology.lyr")
            arcpy.ApplySymbologyFromLayer_management(wclayertemp, wcsymbology)
            arcpy.SaveToLayerFile_management(wclayertemp, wclayer, "RELATIVE", "10.1")



            try:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
                addlayer = arcpy.mapping.Layer(wclayer)
                arcpy.mapping.AddLayer(dataFrame, addlayer, "TOP")
                arcpy.RefreshTOC()
                arcpy.RefreshActiveView()
            except:
                arcpy.AddWarning("Error adding layer to dataframe. Layer can be added from output directory")

                pass

            

            arcpy.Delete_management("in_memory")
            for rasters in inmemorylist:
                arcpy.Delete_management(rasters)


######################### MSAVI Processing End
        #Image processing if WC image already exists in the image archive
        else:
            arcpy.AddMessage(WCimagename + " exists in archive.")

            if ROIshp is None:

                wcmosaic =  "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + ".tif"
                wcmosaicout =  os.path.join(outdir, wcmosaic)
                arcpy.CopyRaster_management(WCarchiveimage, wcmosaicout)

                wclayertemp = os.path.join("WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear )
                wclayer = os.path.join(outdir, "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear )

                metanamearchive = os.path.join(metadir, "meta_" + mlraID + "_" + str(yearID) + ".txt")
                metaname = os.path.join(outdir, "meta_" + mlraID + "_" + str(yearID) + ".txt")
                shutil.copy(metanamearchive, metaname)




            # subset image if region of interest exists
            else:
                landsat = arcpy.sa.Raster(WCarchiveimage)
                spatial_ref = arcpy.Describe(landsat).spatialReference
                arcpy.env.outputCoordinateSystem = spatial_ref

                ROIname = os.path.splitext(os.path.basename(ROIshp))[0]
                wcmosaic = "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "_" + ROIname + ".tif"
                wcmosaicout = os.path.join(outdir, wcmosaic)

                wclayertemp = os.path.join("WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "_" + ROIname )
                wclayer = os.path.join(outdir, "WC_MLRA_" + mlraID + "_" + startyear + "_" + endyear + "_" + ROIname  )

                metanamearchive = os.path.join(metadir, "meta_" + mlraID + "_" + str(yearID) +  ".txt")
                metaname = os.path.join(outdir, "meta_" + mlraID + "_" + str(yearID) + "_" + ROIname + ".txt")
                shutil.copy(metanamearchive, metaname)
                with open(metaname,"a") as metadatac:
                    metadatac.write("\n")
                    metadatac.write("*******************************************************************************" )
                    metadatac.write("\n")
                    metadatac.write("Image was subset from an archived image using: " + ROIname)
                    metadatac.seek(0,0)
                    metadatac.close

                # Check for intersection of inputs
                
                arcpy.AddMessage("Subsetting region of interest from shapefile: " + "\n" + str(ROIshp) + "\n")

                inttest = os.path.join("in_memory", "inttest")
                arcpy.Intersect_analysis([mlradata, ROIshp], inttest)
                rowcount = int(arcpy.GetCount_management(inttest).getOutput(0))

                if rowcount == 0:
                    arcpy.AddError("Error subsetting MLRA " + mlraID + "." + "\n" + "Ensure that the area of interest shapefile intersects the MLRA.")
                    sys.exit(0)
                else:


                    mlradataROI = os.path.join("in_memory", "MLRA_Data_ROI")
                    arcpy.Clip_analysis(mlradata, ROIshp, mlradataROI)
                    cliplayerROI = os.path.join("in_memory","Clip_Layer")


                    arcpy.MakeFeatureLayer_management(mlradataROI, cliplayerROI)

                    clipenvelopeROI = os.path.join("in_memory", "Clip_Envelope")
                    arcpy.FeatureEnvelopeToPolygon_management(cliplayerROI, clipenvelopeROI, "SINGLEPART")
                    descclipROI = arcpy.Describe(clipenvelopeROI)
                    extentclipROI = str(descclipROI.extent).translate(None, 'NaN')
                    

                    arcpy.Clip_management(WCarchiveimage, extentclipROI, wcmosaicout, cliplayerROI,  "255", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
            # Check for no data images
            try:
                MaxWCResult = arcpy.GetRasterProperties_management(wcmosaicout, "MAXIMUM")
                MaxWC = MaxWCResult.getOutput(0)
                arcpy.AddMessage("Max WC is: " + MaxWC)
                DataExists = True
            except Exception:
                DataExists = False
                arcpy.AddWarning("Subset area contains no data.")
                arcpy.Delete_management(wcmosaicout)
                arcpy.AddWarning("No image was produced.")
                pass


            # create layer and symbology
            if DataExists == True:
                arcpy.MakeRasterLayer_management(wcmosaicout, wclayertemp, "#", "#", "1")
                wcsymbology = os.path.join(tooldatapath, "wcsymbology.lyr")
                arcpy.ApplySymbologyFromLayer_management(wclayertemp, wcsymbology)
                arcpy.SaveToLayerFile_management(wclayertemp, wclayer, "RELATIVE", "10.1")

                # Add data to data frame
                try:
                    mxd = arcpy.mapping.MapDocument("CURRENT")
                    dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
                    addlayer = arcpy.mapping.Layer(wclayer + ".lyr")
                    arcpy.mapping.AddLayer(dataFrame, addlayer, "TOP")
                    arcpy.RefreshTOC()
                    arcpy.RefreshActiveView()
                except:
                    arcpy.AddWarning("Error adding layer to dataframe. Layer can be added from output directory")
                    pass

        ### Generate Meta data polygon
        lines = []
        with open(metaname) as file:
            for line in file: 
                line = line.strip() #or some other preprocessing
                lines.append(line) #storing everything in memory!


        Images = []
        WRSID = []
        Date = []
        Year = []


        imageString = 'Total disturbed pixels '

        for line in lines:
            if imageString in line:
                Images.append(line)
                WRSline = line[4:10]
                Yearline = line[10:14]
                Monthline = line[14:16]
                Dayline = line[16:18]
                DateFormated = Monthline + r'/' + Dayline + r'/' + Yearline

                WRSID.append(WRSline)
                Date.append(DateFormated)
                Year.append(Yearline)

        Tiles = list(set(WRSID))
        Tiles.sort()

                
        DatesList = []
        NumImgList = []
        YearsList = []
        NumYearList = []
        for tiles in Tiles:
            DateStringList = []
            YrStringlist = []
            for counter, pathrow in enumerate(WRSID):
                if pathrow == tiles:
                    DateStringList.append(Date[counter])
                    YrStringlist.append(Year[counter])
            
            YrStringlist = (list(set(YrStringlist)))
            YrStringlist.sort()
            YearNum = str(len(YrStringlist))
            ImgNum = str(len(DateStringList))
            DateFullString = "\n".join(DateStringList)
            DatesList.append(DateFullString)
            NumImgList.append(ImgNum)
            NumYearList.append(YearNum)
            Yearfullstring = "\n".join(YrStringlist)
            YearsList.append(Yearfullstring)
        
        



        tempShapeName = os.path.join(outdir, "temp_" + mlraID + "_" + yearID + ".shp")
        arcpy.CopyFeatures_management(mlradata, tempShapeName)
        
        with arcpy.da.UpdateCursor(tempShapeName, "MLRA_ID") as cursor:
            for row in cursor:
                if row[0] != mlraID:
                    cursor.deleteRow()
      
        metaShapeName = os.path.join(outdir, "meta_" + mlraID + "_" + yearID + ".shp")
        arcpy.Sort_management(tempShapeName, metaShapeName, [["Path_Row", "ASCENDING"]])
        arcpy.Delete_management(tempShapeName)

        FieldNames = [f.name for f in arcpy.ListFields(metaShapeName)]
        DelFields = FieldNames[4:]
        arcpy.AddMessage(" , ".join(DelFields))

        arcpy.AddField_management(metaShapeName, "Img_Dates", "TEXT", "", "", "50", "", "NULLABLE", "")
        arcpy.AddField_management(metaShapeName, "Num_Dates", "TEXT", "", "", "50", "", "NULLABLE", "")
        arcpy.AddField_management(metaShapeName, "Img_Years", "TEXT", "", "", "50", "", "NULLABLE", "")
        arcpy.AddField_management(metaShapeName, "Num_Years", "TEXT", "", "", "50", "", "NULLABLE", "")

        
        arcpy.DeleteField_management(metaShapeName, DelFields)
        
        WRS_Field = []
        #Tiles and DatesList
        cursor = arcpy.SearchCursor(metaShapeName)
        for row in cursor:
            WRS_Field.append(row.getValue("Path_Row"))
            arcpy.AddMessage(row)
        arcpy.AddMessage(" , ".join(WRS_Field))
        
        ct = 0
        with arcpy.da.UpdateCursor(metaShapeName, "Img_Dates") as cursor:
            for counter, row in enumerate(cursor):
                
                try:
                    if WRS_Field[counter] in Tiles:
                        row[0] = DatesList[ct]
                        cursor.updateRow(row)
                        ct = ct + 1
                    else:
                        row[0] = "No Data"
                        cursor.updateRow(row)
                except:
                    row[0] = "No Data"
                    cursor.updateRow(row)

        ct = 0
        with arcpy.da.UpdateCursor(metaShapeName, "Num_Dates") as cursor:
            for counter, row in enumerate(cursor):
                try:
                    if WRS_Field[counter] in Tiles:
                        row[0] = NumImgList[ct]
                        cursor.updateRow(row)
                        ct = ct + 1
                    else:
                        row[0] = "No Data"
                        cursor.updateRow(row)
                except:
                    row[0] = "No Data"
                    cursor.updateRow(row)

        ct = 0
        with arcpy.da.UpdateCursor(metaShapeName, "Img_Years") as cursor:
            for counter, row in enumerate(cursor):
                try:
                    if WRS_Field[counter] in Tiles:
                        row[0] = YearsList[ct]
                        cursor.updateRow(row)
                        ct = ct + 1
                    else:
                        row[0] = "No Data"
                        cursor.updateRow(row)
                except:
                    row[0] = "No Data"
                    cursor.updateRow(row)

        ct = 0
        with arcpy.da.UpdateCursor(metaShapeName, "Num_Years") as cursor:
            for counter, row in enumerate(cursor):
                
                try:
                    if WRS_Field[counter] in Tiles:
                        row[0] = NumYearList[ct]
                        cursor.updateRow(row)
                        ct = ct + 1
                    else:
                        row[0] = "No Data"
                        cursor.updateRow(row)
                except:
                    row[0] = "No Data"
                    cursor.updateRow(row)

            metalayertemp = os.path.join("meta_MLRA_" + mlraID + "_" + startyear + "_" + endyear )
            metalayer = os.path.join(outdir, "meta_MLRA_" + mlraID + "_" + startyear + "_" + endyear )
        # create layer and symbology
            
            arcpy.MakeFeatureLayer_management(metaShapeName, metalayertemp)
            metasymbology = os.path.join(tooldatapath, "MLRA_MetaData_Layer.lyr")
            arcpy.ApplySymbologyFromLayer_management(metalayertemp, metasymbology)
            arcpy.SaveToLayerFile_management(metalayertemp, metalayer, "RELATIVE", "10.1")

                # Add data to data frame
            try:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
                newlayer = arcpy.mapping.Layer(metalayer + ".lyr")
                arcpy.mapping.AddLayer(dataFrame, newlayer, "TOP")
                    
                layer = arcpy.mapping.ListLayers(mxd, "")[0] 
                arcpy.AddMessage(layer)
                if layer.supports("LABELCLASSES"):
                    for lblclass in layer.labelClasses:
                            #lblClass.expression = lblClass.expression = "\"<FNT size = '18'>\" & [Img_Years] & \"</FNT>\""
                        lblclass.showClassLabels = True
                        lblclass.expression = "[Img_Years]"
                    
                    #lblclass.expression = '"%s" & [Img_Years] & "%s"' % ("<FNT size=""18"">", "</FNT>")
                        lblclass.expression = '"{}" + [Img_Years] +  "{}"'.format("<FNT name='Arial' size = '12'>""<BOL>","</BOL>""</FNT>") 
                    
                        layer.showLabels = True
                    
                        #Adding Labels layer = arcpy.mapping.ListLayers(mxd, "")[0]  if layer.supports("LABELCLASSES"):     for lblclass in layer.labelClasses:         lblclass.className = "Flaeche"         lblclass.expression = ""<CLR red='255' green='0' blue='0'>" & [FID] & VBNewLine & Round([AREA_mm2],2) & "</CLR>""     lblclass.showClassLabels = True layer.showLabels = True arcpy.RefreshActiveView()  

                arcpy.RefreshActiveView()

                    
            except:
                arcpy.AddWarning("Error adding layer to dataframe. Layer can be added from output directory")
                pass    
        PackageOut = wclayer + '.lpk'
        PackageLayers = [addlayer, layer] 
        arcpy.PackageLayer_management(PackageLayers, PackageOut, "PRESERVE", "CONVERT_ARCSDE", "#", "ALL", "ALL", "CURRENT", "#","RaBETv7","RaBETv7")           


        #arcpy.CheckInExtension("Spatial")


        arcpy.AddMessage("\n" + "Woody cover raster creation is complete for MLRA " + mlraID + " - composite year " + str(endyear) + "." + "\n" )
        return

####################################################### Next Tool



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
        param4.value = os.path.join(sys.path[0], "tooldata","UI.shp")
        #param3.filter.type = "ValueList"
        #param3.filter.list = []

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
        import arcpy


        if str(parameters[3].value) == "Input Shapefile":
            parameters[4].enabled = "false"
            parameters[5].enabled = "true"
            parameters[6].enabled = "true"
            #parameters[4].value = "null"

            if parameters[5].altered:
                #field_names = []

                #parameters[6].filter.list = field_names
                featureclass = parameters[5].value
                #parameters[8].value = str(featureclass)
                field_names = sorted([f.name for f in arcpy.ListFields(featureclass)])
                #field_names.insert(0, "No Name Field")
                parameters[6].filter.list = field_names
        elif (parameters[3].value) == "Draw Polygon on Map Interactively":
            parameters[4].enabled = "true"
            parameters[5].enabled = "false"
            parameters[6].enabled = "false"
            #parameters[4].value = os.path.join(sys.path[0], "tooldata","UI.lyr")
            parameters[5].value = os.path.join(sys.path[0], "tooldata","UI.shp")







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

        # Import modules
        import os
        import time
        import sys
        import re, locale
        import arcpy
        import glob
        import webbrowser
        from operator import itemgetter
        from arcpy import env
        from arcpy.sa import *

        # Define input parameters
        rasterstring = parameters[0].valueAsText.replace("'","")
        rasterlist = sorted(rasterstring.split(";")) # Create and sort list of input rasters
        outdirmain = parameters[1].valueAsText # Output directory - all files output here
        outname = parameters[2].valueAsText # rootname for output files
        method = parameters[3].valueAsText #
        featureset = parameters[4].valueAsText #
        zoneshp = parameters[5].valueAsText # Shapefile input
        zonename = parameters[6].valueAsText # Field name containg unique display name

        arcpy.AddMessage(method)
        charttype = parameters[7].valueAsText # Chart output
        excelfile = parameters[8].valueAsText # excel sheet output
        nodata = parameters[9].valueAsText # NoData handling in zonal statistics

        outdirname = "Analysis_" + outname
        outdir = os.path.join(outdirmain, outdirname)
        if not arcpy.Exists(outdir): 
            arcpy.CreateFolder_management(outdirmain, outdirname)
        else:
            arcpy.AddError("Output directory exists. Please enter a different output name.")
            sys.exit(0)

        if nodata == "true":
            datatype = "DATA"

        else:
            datatype = "NODATA"

        # Check out the Spatial Analyst extension
        arcpy.CheckOutExtension("Spatial")

        #Determine directories containg Python toolbox and reource data
        rootpath = sys.path[0]
        scratchws = os.path.join(rootpath, r"scratch\scratch.gdb")
        tooldatapath = os.path.join(rootpath, r"ToolData")
        mlradata = os.path.join(tooldatapath, r"MLRA_Data.shp")

        # Set Geoprocessing environments
        #arcpy.env.scratchWorkspace = scratchws
        arcpy.env.overwriteOutput = True
        #arcpy.env.workspace = outdir
        landsat = arcpy.sa.Raster(rasterlist[0])
        spatial_ref = arcpy.Describe(landsat).spatialReference
        arcpy.env.outputCoordinateSystem = spatial_ref

        # Output names
        outtablename = outname + ".dbf"
        outtable = os.path.join(outdir, outtablename) # Table output name
        outshapename = outname + "_polygons"
        outshape = os.path.join(outdir, outshapename + ".shp")
        outlayer = os.path.join(outdir, outshapename + ".lyr")
        zonelayer = os.path.join("in_memory", outname + "_zonelayer")
        zonethresh = 10 # max number of graphs that can dispayed before visualization is disabled.

        #Copy input shapefile to outputdirectory and create layer in memory
        if method == "Input Shapefile":
            arcpy.CopyFeatures_management(zoneshp, outshape)
            arcpy.arcpy.MakeFeatureLayer_management(outshape, zonelayer)
            

            if  zonename != None:
                zonefield = zonename
                arcpy.AddMessage("UI fieldname: " + zonefield)
                TitleString = ""
            else:
                fieldnames = [field.name for field in arcpy.ListFields(zonelayer)]
                zonefield = fieldnames[0] # Set the Zone Field to the first OID field
                arcpy.AddMessage("UI fieldname: " + zonefield)
                TitleString = "Polygon "

            #arcpy.arcpy.MakeFeatureLayer_management(outshp, zonelayer)


        elif method == "Draw Polygon on Map Interactively":
            arcpy.AddMessage("feature set")

            arcpy.CopyFeatures_management(featureset,outshape)
            arcpy.arcpy.MakeFeatureLayer_management(outshape, zonelayer)
            arcpy.SaveToLayerFile_management(zonelayer, outlayer)

            fieldnames = [field.name for field in arcpy.ListFields(zonelayer)]
            zonefield = fieldnames[0] # Set the Zone Field to the first OID field
            arcpy.AddMessage("UI fieldname: " + zonefield)

            TitleString = "Polygon "
        else:
            arcpy.AddMessage("Area of Interest Method Not Supported")



        # Add Zone ID field
        arcpy.AddField_management(zonelayer, "Zone_ID", "TEXT")
        expression = r"!" + zonefield + r"!"
        arcpy.CalculateField_management(zonelayer, "Zone_ID", expression, "PYTHON_9.3")

        #Get geometry of MLRA for intersection check
        firstraster = str(os.path.basename(rasterlist[0]))
        mlraID= firstraster.split("_")[2]
        arcpy.AddMessage(mlraID)
        expression = """"MLRA_ID" = '{0}'""".format(mlraID)
        cliplayertile = os.path.join("in_memory", "Clip_Layer_" + outname)
        arcpy.MakeFeatureLayer_management(mlradata, cliplayertile, expression)

        #Check for more than one MLRA in WC rasters
        MLRAList = []
        for rasters in rasterlist:
            MLRAList.append(os.path.basename(rasters).split("_")[2])
        UniqueMLRAs = set(MLRAList)
        arcpy.AddMessage(UniqueMLRAs)
        if len(UniqueMLRAs) > 1:
            arcpy.AddError("WC rasters from multiple MLRAs were input. All input WC rasters must be from the same MLRA.")
            sys.exit(0)

        #Check for intersection of MLRA  boundary and area of interest polygons
        inttest = os.path.join("in_memory", "inttest")
        arcpy.Intersect_analysis([cliplayertile, zonelayer], inttest)
        polyrows = int(arcpy.GetCount_management(zonelayer).getOutput(0))
        rowcount = int(arcpy.GetCount_management(inttest).getOutput(0))
        if rowcount < polyrows:
            arcpy.AddError("Error subsetting MLRA " + mlraID + "." + "\n" + "Ensure that the area of interest polygons are located inside the MLRA.")
            sys.exit(0)

        # Add RaBET_ID field
        fields = arcpy.ListFields(zonelayer)
        fieldnames = []
        for names in fields:

            fieldnames.append(str(names.name))
        if "RaBET_ID" not in fieldnames:
            arcpy.AddField_management(zonelayer, "RaBET_ID", "TEXT")

        if "RaBET_NUM" not in fieldnames:
            arcpy.AddField_management(zonelayer, "RaBET_NUM", "LONG")

        if "RaBET_NAME" not in fieldnames:
            arcpy.AddField_management(zonelayer, "RaBET_NAME", "TEXT")

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
            arcpy.SelectLayerByAttribute_management(zonelayer, "NEW_SELECTION", expression)

            namecursor = str([row[0] for row in arcpy.da.SearchCursor(zonelayer, "RaBET_NAME")][0])

            arcpy.AddMessage("\n" + "Processing statistics for: {0}".format(str(namecursor)))

            for index, rasters in enumerate(rasterlist):
                rastername = os.path.splitext(os.path.basename(rasters))[0]
                Yearname = rastername[-4:]

                currenttable = os.path.join("in_memory" ,str(zones) + Yearname[-2:]   )
                
                
                try:
                    
                    ZonalStatisticsAsTable(zonelayer, "RaBET_ID", rasters, currenttable, datatype, "ALL")
                    arcpy.AddField_management(currenttable, "CV", "DOUBLE")
                    #arcpy.AddField_management(currenttable, "Difference", "DOUBLE")
                    arcpy.AddField_management(currenttable, "Image", "TEXT")
                    arcpy.AddField_management(currenttable, "Year", "TEXT")
                    arcpy.AddField_management(currenttable, "RaBET_Name", "TEXT")

                    arcpy.DeleteField_management(currenttable, "ZONE-CODE")

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

            rowcount = str(int(arcpy.GetCount_management(tables).getOutput(0)))

            if rowcount != '1':
                removetables.append(tables)
            else:
                fieldMappings.addTable(tables)

        #Exclude tables with no data in merge
        for nulltables in removetables:
            tablelist.remove(nulltables)

        #Merge zonal statistics tables in to final output table
        allfieldslist = []
        for field in fieldMappings.fields:
            allfieldslist.append(field.name)

        if field.name not in allfieldslist:
            fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex(field.name))

        arcpy.Merge_management(tablelist, outtable, fieldMappings)

        meanvalues = []
        differencevalues = []
        titlevalues = []


        templayer = outname + "_polygons"
        arcpy.MakeFeatureLayer_management(outshape, templayer)
        arcpy.SaveToLayerFile_management(templayer, outlayer, "RELATIVE")
        arcpy.Delete_management(zonelayer)


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

                    #row.setValue("Difference", float(difference))

                    count += 1
                cursor.updateRow(row)
                meanvalues.append(mean)
                differencevalues.append(abs(difference))
                del row
            del cursor

        #Output Excel spreadsheet
        if excelfile == "true":
            try:
                arcpy.AddMessage("\n" + "Exporting Excel spreadsheet." + "\n")
                excelout = os.path.join(outdir,outname + ".xls")
                arcpy.TableToExcel_conversion(outtable, excelout)
            except:
                arcpy.AddError("Error exporting Excel spreadsheet.")
                pass

        # Output Chart
        if charttype == "true":

            try:
                #arcpy.env.addOutputsToMap = False
                arcpy.AddMessage("\n" + "Creating charts for output." + "\n")
                for zones in uniquezones:

                    expression = """"RaBET_ID" = '{0}'""".format(zones)
                    tview = os.path.join("in_memory", arcpy.CreateUniqueName("tv"))
                    arcpy.MakeTableView_management(outtable, tview, expression)
                    expression = """"RaBET_ID" = '{0}'""".format(zones)
                    #yearcount = str(int(arcpy.GetCount_management(tview).getOutput(0)))
                    try:
                        title = str([row[0] for row in arcpy.da.SearchCursor(tview, "RaBET_NAME", expression)][0])
                    
                        arcpy.AddMessage("\n" + "Creating chart for: " + title)
                        #arcpy.AddMessage("yearcount")
                        #arcpy.AddMessage(yearcount)

                        graphname = outname + "_Chart_" + zones
                        graphimage = os.path.join(outdir, graphname + ".jpg")
                        maxmean = max(meanvalues)

                        if (maxmean >= 0) and (maxmean < 25):
                            graphtemplate = os.path.join(tooldatapath,"GraphTemplateMean25.tee")

                        elif (maxmean >= 25) and (maxmean < 50):
                            graphtemplate = os.path.join(tooldatapath,"GraphTemplateMean50.tee")

                        elif (maxmean >= 50) and (maxmean <= 100):
                            graphtemplate = os.path.join(tooldatapath,"GraphTemplateMean100.tee")

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
                        arcpy.MakeGraph_management(graphtemplate, graph, graphname)

                        # Save the graph as an image
                        if not os.path.isfile(graphimage):
                            arcpy.SaveGraph_management(graphname, graphimage, "MAINTAIN_ASPECT_RATIO", "600", "375")

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
            arcpy.AddMessage("\n" + "Number of polygons exceeds threshold of " + str(zonethresh) + "." + "\n" + "Charts will not be displayed, but were saved to output directory." +"\n" + "\n")

        arcpy.Delete_management("in_memory")

        dropFields = ["Zone_ID", "RaBET_ID", "RaBET_NUM"]
        arcpy.DeleteField_management(outlayer, dropFields)

        # Add polgons to the data frame
        try:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
            addlayer = arcpy.mapping.Layer(outlayer)
            arcpy.mapping.AddLayer(dataFrame, addlayer)
            #addlayer.labelClasses[0].expression = "[RaBET_NAME]"
            #addlayer.showLabels = True
            arcpy.RefreshActiveView()
            arcpy.RefreshTOC()
        except:
            arcpy.AddWarning("Error adding layer to dataframe. Layer can be added from output directory")
            pass

        return
