NPRE 555 CP1 “Monte Carlo Code of Neutron Transport in a 2D slab”
By Muhammad Abdelghany
PhD Student at the NPRE department at the UIUC
***********************************************

This code is made using MATLAB and is composed of two main files.
The first file is called “cp1_cyc.m”. This is the main file, which is used to make a run
for a specific number of cycles that is required to be performed to carry out our calculation
of the multiplication factor (Keffective), and the flux distribution over the width of the
slab.

The second file is a function called “ke” which is the main core of the code. This file
returns the value of the multiplication factor for each cycle, as we run the main file
 “cp1_cyc.m”.

The Details of Each File:
*************************

1. “cp1_cyc.m”
--------------

This file starts with specifying the number of cycles needed for our run, then it loops 
over the function “ke” which return the value of “k_effective” for each cycle, as 
calculated by the second file. Each vale of “k_effective” is stored in a row matrix 
called “m” The size of this loop is equal to the inputed number of cycles. Then it 
calculate the average value of K_effective over all the cycles by taking the average 
of the elements in the matrix “m”.

2. “ke.m”
---------

This file starts with specifying the number of particles through the variable “npart”,
which can be changed according to the needed accuracy. Next, it starts a for-loop
on theses particles, inside each of them we track the path of the neutrons inside the
slab. For each particle “j”, its position is started at x(j), y(j) by generating 
a random variable that is uniformly distributed between [0Cm , 100Cm] for x(i) and between 
[0Cm , 1Cm] for y(i). “x” and “y” are treated differently because “x” is finite, but “y”
is infinite.
After that it checks whether the particle is born at the left region or the right region,
and specify which cross-sections it should use.
Next, it samples randomly the travele distance “r”, and it also generates a 
uniform random number for the angle “theta” between [0, 2*pi], and then calculate for the
new position. After that it checks whether the particle is still in the first or the second
region, and according to that check wether the particle will be absorbed or not. If
absorbed it will generate k_infinity neutrons, which then is assigned as the flux 
contribution of particle “j” at the position x(j) , y(j). If the particle neither absorbed
nor leaked form the system, the loop will continue to track the new positions of the 
particle due to the scattering interaction, till it is absorbed.
the code also have a counter “history(j)” that calculate the number of scattering events
happend before the particle is absorbed. At the end, it calculates the flux at each position
“x” due to the contribution from all the particles. At this point we should discretize the space 
“x” into intervlas, I choose them to be 100 intervals after making some iterations to check
the a convienient, in order to fit the calculated flux at each space mesh.
We draw the flux distribution from all particles as a function of “x”.
We also calculate the integrated fux “alpha” all over “x” and find the total number of 
neutrons exist in the slab, which can be used to calculate the multiplication factor, by 
dividing this integrated flux by the initial number of particle “npart”.


The Additional Folder “CP1 without Cycles”
******************************************

It include two codes, exactly similar to this code but with removing the option of 
carrying out the calculations for multiple cycles, but it is only for one cycle.
One of these codes have the slab is discretized into 100 meshes, and the other is 
discretized into 1000 meshes, while calculating the flux distribution. This is to study
the effect of refining the discretization of the space and check whether this will enhance
the shape of the final flux distribution or not.
