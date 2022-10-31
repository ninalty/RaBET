from RaBET.wcmaptool import WCMapTool

kk = WCMapTool()

imagedir = 'C:\\Users\\ninal\\Documents\\postdoc UCDavis\\testt\\data_Chandra\\MLRA 41_test\\LC080360372015061701T1-SC20210721135904'  # Landsat image directory
mlraID = '41'  # MLRA name/ID
yearID = str(2015)  # Year to process
outdir = 'C:\\Users\\ninal\\Documents\\postdoc UCDavis\\testt\\output'  # Qutput directory
ROIshp = None #'C:\\Users\\ninal\\Documents\\postdoc UCDavis\\testt\\data_Chandra\\tool_data\\MLRA_Data.shp'  # Optional ROI shapefile
parameters = [imagedir, mlraID, yearID, outdir, ROIshp]

kk.execute(parameters)