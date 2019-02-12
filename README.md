# Codes for the STIMULUS paper
- GenerateFactorTable.m pre-calculates a look-up table of AFDs given electrode locations; run this first before STIMULUS_all*.m
- STIMULUS_all*.m calculates results using STIMULUS
- TI_resistors.m caculates TI stimulation results
- generage_figure.m loads pre-calculated results and generates Fig.8-Fig.12
- Auxiliary functions that implement parts of the algorithm
- HodgkinHuxleyFast.cpp: mex file for classic HH equations, using the Boost library. Precompiled on Linux
- HodgkinHuxleyE_iso.cpp: mex file for excitatory HH equations, using the Boost library. Precompiled on Linux
