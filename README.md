[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Python-version:](https://img.shields.io/badge/python-v3.7-green)]() | [FitViewer](#fitviewer) | [Usage](#usage) | [Requirements](#requirements) | [References](#references) 

# FitViewer
[FitViewer](https://github.com/elija-feigl/FitViewer) for context based zoning and cropping of cryo-EM maps using pseudoatomic models fitting into cryo-EM electrondensity maps of lattice based DNA origami structures

# Usage
start the interactive Jupyter-Notebook by executing the the following command from the repository directory:
```
jupyter notebook viewertool.ipynb
```
Step by step instructions are provided in the Jupyter notebook 

# Requirements
[FitViewer](https://github.com/elija-feigl/FitViewer) is currently distributed as a [Jupyter Notebook](https://jupyter.org/documentation). It requires an [ipytohn kernel](https://ipython.readthedocs.io/en/stable/install/kernel_install.html) with the follwoing dependencies.

dietzlab_Nanodesign;
custom python3 version of [Nanodesign by Autodesk](https://github.com/Autodesk/nanodesign)
available at https://github.com/elija-feigl/nanodesign_dietz

required PyPI packages:
```
ipywidgets >= 7.5.1
mrcfile >= 1.1.2
MDAnalysis >= 1.0.0
PyQt5 >= 5.14.2
nglview >= 2.7.7
numpy >= 1.18.4
attrs >= 19.3.0
```

## Supported formats
[FitViewer](https://github.com/elija-feigl/FitViewer) relies on the [MDAnalysis](http://www.mdanalysis.org/) package for reading trajectories or coordinate files,
[mrcfile](https://github.com/ccpem/mrcfile) for handling cryo-Em data, and DNA Origami design files generated using [cadnano](https://cadnano.org)

# References
When using [FitViewer](https://github.com/elija-feigl/FitViewer) in published work, please cite the following paper:

manuscript status:  10.07.2020 currently under review

