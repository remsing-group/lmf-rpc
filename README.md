# lmf-rpc

LAMMPS pair styles for Gaussian truncated electrostatics, with cutoff Lennard-Jones for van der Waals interactions. 
The pair style is lj/cut/coul/trunc. 
In the LAMMPS input file, use "pair_style	lj/cut/coul/trunc <sigma> <LJ cutoff> <GT cutoff>"
LJ epsilon and sigma parameters are added in pair_coeff commands as usual. 
Note that the pair style has not been modified to work with RESPA. 

This pair style is for LAMMPS July 30, 2016 version. 
To include, add the .cpp and .h files to the src directory and compile LAMMPS as usual. 
