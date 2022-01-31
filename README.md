# 16S-rRNA-pipeline-modules
Collection of modules made during (Graduation) internship .This repository contains the finalized modules as implemented in the pipeline, as well as modules used during testing and validation (including case study). 

Main modules of workflow:
- sankey-file-prep.py
- sankey-diagram-html-generator.R
- sankey-diagram-png-generator.R

"Metataxonomics" folder contains the pipeline as a whole. The newly developed Docker file is stored at "metataxonomics/Docker/sankey-dockerinstall-Rpackages"

"Others" folder contains modules i used for the case study and for testing/development (overlay, grid)

The most important module in the folder "Others" is: screenshot-html-files-Chrome.py   ; this was used for the output (.html files sankey diagram) 
This module uses the .html files in the "html files" folder to make .png files in which the names are to the left of the nodes, adds an overlay and is able to form a grid of various sankey diagrams. 

Things to note:
- Path to folders are still hardcoded, change manually
- Make sure there is a folder to store the .html files (folder: "html files") and .png files (folder: "results - png sankey")
- You have Chrome as a browser and the chromedriver (chromedriver.exe) correspondent to that version of chrome. https://chromedriver.chromium.org/getting-started
