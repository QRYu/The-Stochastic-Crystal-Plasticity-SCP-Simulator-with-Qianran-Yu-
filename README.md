# The-Stochastic-Crystal-Plasticity-SCP-Simulator-with-Qianran-Yu-
The SCP simulation package contains a computer code written in C++. It is a newly developed stochastic solver for crystal plasticity problem. We have confirmed its validity with the case of uniaxial loading of tungsten single crystals. 

****Developed by****

Qianran Yu

It is an attempt to use residence time algorithm to solve tradiational crystal plasticity problem with rate theory method. Please send an email to the following address:

Qianran Yu (yuqianran0709@gmail.com)

Jaime Marian (jmarian@ucla.edu)

****Introduction****

The traditional crystal plasticity theory assumes that slip systems are independent of each other such that they, during straining, are uniformly activated. However, such treatment is at odds with physical observation. In this model, we regarded shear rates as "reaction rates", and used the general rate theory method to solve the CP problem. The stochasticity of slip system motions are realized by residence time algorithm. Please see the following journal paper for detail:

[1] Qianran Yu, Enrique Martinez, Javier Segurado, Jaime Marian, “A stochastic solver based on the residence time algorithm for crystal placticity models,” Computational Mechanics 68, 1369-1384 (2021). (link: https://link.springer.com/article/10.1007/s00466-021-02073-7)

****How to use****

Before use:

Please install gcc/g++ for the newest version (c++11 or newer).
Please either unzip eigen3 library and put it into src directory, or install eigen3 in local computer.

Open terminal and type "make" under "src" directory, then and executable file named "cpexe" will be generated. Run the simulations using command "./cpexe".


