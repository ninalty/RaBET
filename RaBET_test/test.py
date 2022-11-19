from RaBET.wcmaptool import WCMapTool
import os

kk = WCMapTool()

local_path = 'C:\\Users\\ninal\\Documents\\postdoc_UCDavis\\RaBET_v1\\'

imagedir = os.path.join(local_path, 'RaBET_test\\image_input\\LC08_L1TP_034038_20210619_20210629_02_T1') # Landsat collection2
# imagedir = os.path.join(local_path, 'RaBET_test\\image_input\\LC08_L1TP_034038_20180611_20200831_02_T1') # Landsat collection2
# imagedir = os.path.join(local_path, 'RaBET_test\\image_input\\LC08_L1TP_034038_20150603_20210721_01_T1') # Landsat collection 1

mlraID = '41'  # MLRA name/ID
yearID = str(2021)  # Year to process
outdir = os.path.join(local_path, 'RaBET_test\\output')  # Qutput directory
ROIshp = None # os.path.join(local_path, 'RaBET_test\\tooldata\\ROI.shp') # Optional ROI shapefile
parameters = [imagedir, mlraID, yearID, outdir, ROIshp]

# check whether scratch, output and WC_Image_Archive folders exist.
if not os.path.exists(outdir):
    os.mkdir(outdir)

scratch_dir = os.path.join(local_path, 'RaBET_test\\scratch')
if not os.path.exists(scratch_dir):
    os.mkdir(scratch_dir)

archive_dir = os.path.join(local_path, 'RaBET_test\\WC_Image_Archive')
if not os.path.exists(archive_dir):
    os.mkdir(archive_dir)

kk.execute(parameters, 2)