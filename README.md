# lmf-rpc

There are two LAMMPS pair styles for Gaussian truncated electrostatics, with cutoff Lennard-Jones for van der Waals interactions. 
The first pair style is lj/cut/coul/trunc. 
In the LAMMPS input file, use "pair_style	lj/cut/coul/trunc <sigma> <LJ cutoff> <GT cutoff>"
LJ epsilon and sigma parameters are added in pair_coeff commands as usual. 
Note that the pair style has not been modified to work with RESPA. 

The second pair style is lj/cut/tip4p/trunc and is specific to TIP4P water models. 
In the LAMMPS input file, use "pair_style	lj/cut/coul/trunc <sigma> <oxygen atom type> <hydrogen atom type> <O-H bond type> <H-O-H angle type> <qdist> <ljcutoff> <coulomb cutoff>"
qdist is the distance of the M-site from the oxygen site, as discussed in the LAMMPS manual. 
LJ epsilon and sigma parameters are added in pair_coeff commands as usual. 

Also included is a LAMMPS "fix" for performing simulations in the presence of the LMF potential. 
In the LAMMPS input file, use "fix <name> <group to apply field to> "lmf" <number of points in LMF potential file> <LMF potentail file name>"
The required input is a file containing z-coordinate in column 1, LMF potential in column 2, and force from the LMF potential in column 3. 
The energy units must be consistent with those specified in the LAMMPS input file. 

This pair style is for LAMMPS July 30, 2016 version. 
To include, add the .cpp and .h files to the src directory and compile LAMMPS as usual. 
