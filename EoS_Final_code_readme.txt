Copyright (c) 2018.
This file was created by Paolo Parotto and Debora Mroczek, Department of Physics, University of Houston, Houston, TX 77204, US.


This is a beta version of the BEST Collaboration program producing an EoS matching Lattice QCD at muB=0, and containing a critical point in the 3D Ising universality class, in a parametrized form.
It allows for different choices of constraints on the shape of the critical line, which reduce the number of parameters.


THE INPUTS
The inputs required to run the program are:
- Input Lattice data. They correspond to the 0th, 2nd and 4th coefficients in the
Taylor expansion of the pressure. They are contained in the file "Lattice_Data_5_821_dT1.dat"
- Input HRG pressure data. This is a table containing the values of the pressure, calculated making use of a non-interacting Hadron Resonance Gas (HRG) model in the computational range. The data are contained in the file "Press_HRG_MUB000601_T005300_dT1.dat".
- Input parameter file. This file contains information on how the choice of parameters in the Ising -> QCD map is made. More on this will be in the following.


THE BODY OF THE PROGRAM
NOTE: The computational grid of the program is T = 5 - 821 MeV,
muB = 0 - 601 MeV, with grid size of 1 MeV. The output of the program is in the 
range T = 30 - 800 MeV, muB = 0 - 450 MeV.
The program goes through the following steps:
1.) Import Lattice data from the location stated above, multiply by T^4 and store. 
NOTE: throughout the program, all quantities are dimensional. At the end, thermodynamic quantities exported are normalized with powers of the temperature.
2.) Set the parameters, importing from the parameter file (see below).
3.) Given the choice of parameters, a coordinates file is created. It contains the map between (T,muB) and (R,theta). Before creating the file, the code checks for its presence in the working directory; if it is already present, it imports it and moves on.
4.) Using the coordinates file, the Ising contribution to the expansion coefficients is calculated.
5.) Using the coordinate file, the 3D pressure and its derivatives are calculated at all points in the grid.
NOTE: the calculation makes use of coordinate points with negative muB too, due to the symmetrization procedure described in the paper.
6.) The Non-Ising contribution to the expansion coefficients is calculated from the Lattice coefficients and the Ising ones.
7.) The Non-Ising contribution to the pressure is calculated in all points in the grid, as a Taylor expansion of the Non-Ising coefficients.
8.) The regular pressure is filtered through a Gaussian filter. This is done in order to reduce unphysical wiggles in the thermodynamic functions at largu muB due to the truncation of the Taylor series. The parameters of the filter do not need to be entered manually.
9.) Numerical derivatives of the Non-Ising pressure are calculated.
10.) Ising and Non-Ising contributions to the pressure ad its derivatives are combined (summed).
11.) HRG pressure is imported and multiplied by T^4. Numerical derivatives of the pressure are calculated.
12.) HRG and Ising+Non-Ising contribution to the pressure and derivatives are combined via a smooth merging. The merging makes use of an hyperbolic tangent. The parameters of the merging do not need to be input manually.
13.) Once derivatives are defined, they are combined and normalized to generate the thermodynamics. Pressure, entropy density, baryon density, energy density, speed of sound and second baryon number cumulant are exported in the range T = 30 - 800 MeV, muB = 0 - 450 MeV. They are all normalized with powers of the temperature.


THE PARAMETER FILE
The parameter file needs to contain the "mode" of input, corresponding to the way the location on the critical point is constrained, and then depending on the mode the additional needed parameters. 
1.)
"Mode": FREE
This means the choice of the six parameters in the map is completely free. In this case, the additional parameters needed are: TC, the temperature at the critical point; muBC, the chemical potential at the critical point; angle1, the angle between the h=0 axis and the T=0 axis (in degrees); angle2, the angle between the r=0 axis and the T=0 axis (in degrees); w, the global scaling parameter in the mapping; rho, the relative scaling in the mapping. 
Example: FREE 140 350 6 96 1 2
2.)
"Mode": PAR
This means the critical point will lie on a parabola. In this case, the additional parameters needed are: T0, the value of temperature at which the parabolic pseudo-critical line crosses the T axis; kappa, the curvature of the transition line at the T axis; muBC, the chemical potential at the critical point: since both TC and angle1 are determined by this choice, they do not need to be input; anglediff, the difference between angle2 and angle1 (in degrees), this is usually set to be 90 degrees; w, rho as before. NOTE: the values of T0 and kappa currently used are from http://inspirehep.net/record/1385115.
Example: PAR 155 -0.0149 350 90 1 2
1.)
"Mode": STR
This means the critical point will lie on a straight line. In this case, the additional parameters needed are: T0, the value of temperature at which the straight pseudo-critical line crosses the T axis; muBC, the chemical potential at the critical point; angle1, the angle between the h=0 axis and the T=0 axis (in degrees); angle2, the angle between the r=0 axis and the T=0 axis (in degrees); w, the global scaling parameter in the mapping; rho, the relative scaling in the mapping. 
Example: STR 155 350 3 93 1 2


NOTE ON THE PARAMETERS
This project is intended for use in the framework of the BES-II. The program has been tested for "realistic" values of the parameters.
For the location of the critical point, the code is intended to place it in the chemical potential range spanned by the BES-II program (muB < 500 MeV).
The values of the angles should not be multiples of 90 degrees. This would prevent the critical contribution from the Ising model to be reflected in the QCD thermodynamical quantities that are calculated.
The code was tested for values of w = 0.2 - 2, rho = 0.5 -4. Larger values of w should not result in problems; however, the larger w is, the more the critical contribution is "washed away". 

NOTE ON THE FILE
The public version of the code also contains the tables for the different thermodynamic quantities for the set of parameters included in the manuscript, as well as a different choice with w=0.75.

COMPILING
This code uses the standard libraries <stdio.h>, <stdlib.h>, <math.h>, <string.h>, <time.h>. In addition, the libraries sys/types.h, sys/stat.h and direct.h were included in order to use the functions _mkdir() and _chdir(). Please note that different systems might use the functions mkdir() and chdir() (no underscore), and these might be defined in other header files. Please note that the syntax should be changed if necessary where these functions are called (lines 219, 220, 465, 482).
Please, do not hesitate to contact the authors for questions/problems.

Standard: compile with
gcc -o (name of executable) *.c
or
gcc -lm -o (name of executable) *.c

Try also:
gcc -l/usr/include -c *.c
gcc -L/usr/lib/ *.o -lgsl -lgslcblas -lm -o EOS

NOTE: this is a beta version, for problems with compiling contact the authors.


CONTACT 
For problems, debugging, contact Paolo Parotto at pparotto@uh.edu


WHEN USING THIS CODE
When using this code, the following work should be cited:
http://inspirehep.net/record/1672952