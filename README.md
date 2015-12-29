# Ray based FEM solver for the High-Frequency Helmholtz Equation

This a extremely experimental package that implements a discretization of the Helmholtz equation with a non polynomial basis enriched with a geometric optics Ansatz. The basis are probed using a novel adaptive learning approach that relies on solving a low-frequency problem.  

Moreover, all linear systems (so far only the low-frequency) are solved using a modification of the method of polarized traces for FEM. 

# Dependencies

This Package relies on the iFEM package of professor Long Chen. In order to run the tests you need to have the iFEM package in your MATLAB path. 


Jun Fang, Leonardo Zepeda-Núñez, and Hongkai Zhao
University of California, Irvine 2015