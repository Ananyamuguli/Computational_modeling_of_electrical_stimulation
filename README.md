# Computational modeling and analysis of acute electrical brain stimulation paradigms during Secondary Brain Injury (SBI).

The project uses PyDSTool library (https://pydstool.github.io/PyDSTool/FrontPage.html) to work with the ODEs presented in the paper. To properly work with PyDSTool:
  1. python version <= 3.9.13
  2. setuptools <= 65.0.0
The conda environment that provides these is provided in conda_environment_PyDSTool.yml. To set up this envornment, go to the base environment and execute the following command:
conda env create -f conda_environment_PyDSTool.yml
