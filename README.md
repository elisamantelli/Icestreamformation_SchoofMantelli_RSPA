# Icestreamformation_SchoofMantelli_RSPA_2021

This repository contains the code used to produce the 3D ice stream simulations in "Ice Stream Formation", by C. Schoof and E. Mantelli, to appear in Proceedings of the Royal Society of London A in 2021.

The code can be run in 2 different modes: 

1. The routine script_divide_runs.m calculates velocity, temperature, surface elevation and effective pressure for steady ice sheet flow lines with subtemperate sliding. It uses shallow ice mechanics. This code is used to produce output plotted in fig. 1 in the paper and fig. 2 in the supplementary materials.

2. The routine run_divide_int.m calculates velocity, temperature, surface elevation and effective pressure for a three dimesnional ice sheet "strip" starting from an ice divide, perturbed with white noise in the friction coefficient. This code solves the Stokes problem for the transverse (y,z) velocity field, a shallow ice model for the x-velocity, mass conservation, and a shallow heat equation. Sliding is both temperature- and effective pressure-dependent. This routine is used to produce all the output from 3D simulations displayed in the paper and supplementary materials. 

Full output for the simulations plotted in the paper is available upon request, as are plotting scripts.  Please, do not hesitate to get in touch for any questions/clarifications/.. (Elisa Mantelli, mantelli@princeton.edu or elisa.mantelli@gmail.com).
