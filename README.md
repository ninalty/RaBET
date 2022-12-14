# RaBET

This repository houses a Python toolbox of tools for preprocessing inputs and postprocessing outputs of the RaBET (Rangeland Brush Estimation Toolbox). To find out more about RaBET, please visit https://www.tucson.ars.ag.gov/unit/publications/PDFfiles/2271.pdf.

## Getting Started
* Ensure you have [ArcGIS Pro](https://pro.arcgis.com/en/pro-app/latest/get-started/get-started.htm) installed.
* Download the toolbox and place it in an appropriate folder on your machine. Navigate to the folder in Catalog. If you expand all the toolsets, you will see the following:

  ![arcgis_pro platform](https://user-images.githubusercontent.com/35977606/207509333-7bfd379b-7af8-44ca-b301-d0e3683cfc78.JPG)

* In order to use the tool, you will need to have the following inputs available:
   * step 1 Under the working directory eg., RaBET_version30 with the following sub-folders: MLRA_41_test folder to hold the Landsat imagery; tooldata folder for shapefiles provided in this repository; output folder to save all outputs from RaBET toolbox. 
   ![directory](https://user-images.githubusercontent.com/35977606/207519859-03492c92-c560-426b-84d7-ec7bc01d8ffd.JPG)
   * step 2 download Landsat test images and move it to MLRA_41_test folder.[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7430812.svg)](https://doi.org/10.5281/zenodo.7430812)
   * step 3 run Generate WC Maps tool.
 
   ![process](https://user-images.githubusercontent.com/35977606/207519976-7197616e-331a-4cce-8e58-1d6e2c34c7fe.JPG)

## Development
For any future developments, add new python script to scripts folder and update the .pyt file accordingly.
