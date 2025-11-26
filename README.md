📘 GPS LEO Geometry Analysis
This repository contains two Jupyter notebooks and a supporting Python module used to analyze GPS satellite geometry relative to a Low Earth Orbit (LEO) satellite.

🚀 How to Use
Install the required Python libraries (see below).
Ensure gnss_lib.py is available in your Python path or in the same directory as the notebooks.
Open the notebooks in JupyterLab or Jupyter Notebook.
Run all cells from top to bottom.

📂 Notebooks Overview
1. Testing.ipynb (Main Analysis Notebook)
This is the primary notebook containing all calculations, all analysis, and all results.
It performs:
GPS satellite position computation
LEO satellite trajectory computation
GPS–LEO visibility checks
Dilution of Precision (DOP) / geometry evaluation
Full plotting, analysis, and validation
Calls to helper functions implemented in gnss_lib.py
Run this notebook to reproduce the complete GPS–LEO geometry analysis.

3. leo_constelation.ipynb (LEO Constellation Generator)
This notebook is used only to generate LEO constellations.

📦 Required Libraries
Make sure you have the following core packages:
pip install numpy scipy matplotlib pandas jupyter

📁 Required Local Module
gnss_lib.py --> https://github.com/Stanford-NavLab/gnss_lib_py
The notebooks will not run without this module.
