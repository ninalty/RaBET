# RaBET

This repository houses a Python toolbox of tools for preprocessing inputs and postprocessing outputs of the RaBET (Rangeland Brush Estimation Toolbox). To find out more about RaBET, please visit https://www.tucson.ars.ag.gov/unit/publications/PDFfiles/2271.pdf.

## Tools
* RaBET_version_27 is the original code works exclusively for collection 1 level 2 Landsat images using ArcPY 2.7.
* RaBET_version_30 is the modified code works for both collection 1 level 2 and collection 2 level 1 Landsat images using ArcPY 3.0.
* 
## Getting Started
* Ensure you have [ArcGIS Pro](https://pro.arcgis.com/en/pro-app/latest/get-started/get-started.htm) installed.
* Download the toolbox and place it in an appropriate folder on your machine. Navigate to the folder in Catalog. If you expand all the toolsets, you will see the following:

  ![arcgis_pro platform](https://user-images.githubusercontent.com/35977606/207509333-7bfd379b-7af8-44ca-b301-d0e3683cfc78.JPG)

* In order to use the tool, you will need to have the following inputs available:
   * step 1 Under the working directory eg., RaBET_version30, you should have tooldata folder provided in this repository. 
   ![directory](https://user-images.githubusercontent.com/35977606/207519859-03492c92-c560-426b-84d7-ec7bc01d8ffd.JPG)
   * step 2 download Landsat test images and move it to MLRA_41_test folder. The original data is from [USGS](https://earthexplorer.usgs.gov/). A data sample is provided for practice. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7430812.svg)](https://doi.org/10.5281/zenodo.7430812)
   * step 3 run Generate WC Maps tool.
 
   ![process](https://user-images.githubusercontent.com/35977606/207519976-7197616e-331a-4cce-8e58-1d6e2c34c7fe.JPG)

## Future Development
* For any future developments, add new python script to scripts folder and update the .pyt file accordingly.

