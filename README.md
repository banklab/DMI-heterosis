# DMI-heterosis

Repository for the manuscript **The influence of incompatibilities and heterosis on hybrid population genetics**

File structure:

All scripts are contained in the folder `src`

Files description:

`Heterosis_DMIs.nb/Heterosis_DMIs.pdf` contains the Mathematica notebook/pdf containing the model description and analytical solutions to limiting cases and approximations

`dmi_heterosis_functions.py` contains the discretization of the model and the implementation of three loci for linked neutral variation

`generate*.py` uses the previous functions file and includes all scripts to generate datasets for linked neutral variation

`fit_bgs.R` contains the code used for the model and supplementary figures of recurrent hitchhiking

`plot*.r` includes all functions for all other plots of the manuscript

