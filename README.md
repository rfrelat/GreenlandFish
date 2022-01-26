# Dynamics of Greenland fish communities


This repository contains the scripts and dataset to reproduce the results presented in the manuscript **Deep demersal fish communities respond rapidly to warming in a frontal region between Arctic and Atlantic waters** by Margrete Emblemsvag, Karl Michael Werner, Ismael Núñez-Riboni, Romain Frelat, Helle Torp Christensen, Heino O. Fock and Raul Primicerio published in *Global Change Biology*.



The script is commented as much as possible, but it is not intended to be a tutorial. For tutorial about tensor decomposition, please visit: https://rfrelat.github.io/Multivariate2D3D.html. For more information about the material and method, please see the original manuscript.



The analysis is carried out in the script [TensorGreenland.R](https://github.com/rfrelat/GreenlandFish/blob/main/TensorGreenland.R). All steps needed to run the Tensor decomposition on the log-scaled abundance and the clustering of the species based on their spatio-temporal dynamics. All figures shown in the manuscript and in the supplementary materials are reproduced while running the script.



Packages to run these scripts are: ade4, PTAk, raster, maps, mapdata, and RColorBrewer. If some of these packages are not installed, you can install them with the following command line:

```{r}
install.packages(c("ade4", "PTAk", "raster", "maps", "mapdata", "RColorBrewer"))
```

Version: all analyses were run on R4.0.2.



The dataset and additional functions are all included in the Rdata file [TensorGreenland.Rdata](https://github.com/rfrelat/GreenlandFish/raw/main/TensorGreenland.Rdata) which contains:

- logtensor: the tensor of log transformed abundance of 55 species, over 18 years and 7 depth layers
- biogeography: the biogeography of the 55 species
- sss, ice, slat, sst: the correlation maps of Sea surface salinity, Sea ice concentration, Sea level air temperature and Sea surface temperature, respectively.
- two functions: `myHeatmap`to plot heatmaps and `hc_wss` to compute the within sum of square of hierarchical clustering.





[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)