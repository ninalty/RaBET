from RaBET.wcmaptool import WCMapTool
import os

kk = WCMapTool()

local_path = 'C:\\Users\\ninal\\Documents\\postdoc UCDavis\\RaBET_v1\\'

imagedir = os.path.join(local_path, 'RaBET_data\\MLRA 41_test\\LC080360372015061701T1-SC20210721135904')  # Landsat image directory
mlraID = '41'  # MLRA name/ID
yearID = str(2015)  # Year to process
outdir = os.path.join(local_path, 'RaBET_test\\output')  # Qutput directory
ROIshp = None # Optional ROI shapefile
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

kk.execute(parameters)