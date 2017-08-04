# PQIndices_DCGrids
Power Quality indicators for DC grids (networks)

There is no actual standard for defining PQ in DC networks, besides some normatives from Military applications and Telecom Networks. 
Therefore, this was was one of the first attempts in defing such indicators that could play a role in the design of a DC grid (microgrid), inlcuding the stability analyis of the controllers. 

All indices were defined based on non-stationary signal processing analyis for both time and frequency domains.

# Notes:
All definitions and rationale for their choice can be found in the paper "Analytical derivation of PQ indicators compatible with control strategies for DC microgrids", published at IEEE PowerTEch Conference, Manchester, UK 2017.
How to cite this work: [x]	I. Ciornei, L. Hadjidemetriou, M. Albu, M. Sanduleac, E.Kyriakides, “Analytical derivation of PQ indicators compatible with control strategies for DC microgrids,” IEEE PowerTech, Manchester, UK,pp. 1-6, 18-22 June 2017 [DOI:10.1109/PTC.2017.7981179]. link: https://www.researchgate.net/publication/318494768_Analytical_derivation_of_power_quality_indicators_compatible_with_control_strategies_for_DC_microgrids

# HOW TO USE:
The mat files in the current folder are either simulated signals or lab aquired signals from an experimnetal setup presented in the paper. 
The Matlab simulink files that were used to generate both lab tests signal and the simulated signals are part of Project: SimplifiedDCMicrogridModel.

There are multiple main files, one for each case:
- manually set up the duty clycles in the control loop of the MGs such that to overcome the variation of the load for a targeted output of 36V and 48V, respectivelly. There are 2 cases for each voltage: Load1=114Ohms, and Load2=228Ohms;
- we prefered manual setup instead of automatic closed loop control in order to evaluate the impact of these changes in the aquired signal

One may run independent each case, or make another main where thses cases are presented all together.
