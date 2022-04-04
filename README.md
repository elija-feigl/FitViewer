[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Python-version:](https://img.shields.io/badge/python-v3.8-green)]() | [FitViewer](#fitviewer) | [Usage](#usage) | [Requirements](#requirements) | [References](#references)

# FitViewer
[FitViewer](https://github.com/elija-feigl/FitViewer) for context based zoning and cropping of cryo-EM maps using pseudo-atomic models fitting into cryo-EM electron density maps of lattice based DNA origami structures

# Usage

start the interactive Jupyter-Notebook by executing the following command from the repository directory:
```
jupyter notebook viewertool.ipynb
```
Step by step instructions are provided in the Jupyter notebook


# Requirements
[FitViewer](https://github.com/elija-feigl/FitViewer) is currently distributed as a [Jupyter Notebook](https://jupyter.org/documentation). It requires an [ipytohn kernel](https://ipython.readthedocs.io/en/stable/install/kernel_install.html) with the following dependencies:
```
python >= 3.8.0
```

[dnaFit](https://github.com/elija-feigl/DNA_Fit) DNA Origami atomic model package

required PyPI packages:
```
ipywidget
ipython
jupyter
```

## Supported formats
[FitViewer](https://github.com/elija-feigl/FitViewer) relies on
the [MDAnalysis](http://www.mdanalysis.org/) package for reading trajectories or coordinate files,
[mrcfile](https://github.com/ccpem/mrcfile) for handling cryo-EM data,
and DNA Origami design files generated using [cadnano](https://cadnano.org)

# References
When using [FitViewer](https://github.com/elija-feigl/FitViewer) in published work, please cite the following paper:

Kube, M., Kohler, F., Feigl, E. et al. Revealing the structures of megadalton-scale DNA complexes with nucleotide resolution. Nat Commun 11, 6229 (2020). https://doi.org/10.1038/s41467-020-20020-7

