# Attractor dynamics for networks with learning rules inferred from in vivo data (U. Pereira & N. Brunel, Neuron 2018)

This repository contains three directories:

1. Dynamics
2. CapacityFung
3. MaximalCapacity
 

# Dynamics
It contains the code for reproducing Fig 6 A-F. 
By changing the value of the parameter A from 9.5 to 3.5 and setting the input current as:

1. Standard normal iid for a novel stimulus
2. One of the stored patterns for familiar stimulus,

the results on Fig 3 and Fig 4 A are also reproduced.

For generating fig4.pdf, first run `python main.py` on the `numeric/network_simulation` directory. The simulation
will take a couple of hours depending on your machine. It will create two large files `the_dynamics.p' (~800MB) and `the_overlaps' (~1.15GB) which
store the dynamics and the overlaps corresponding to   10 realizations of the dynamics.

# Capacity_vs_g
It contains the code for reproducing Fig 4 B-C. 
The same code can be used to perform the parameter exploration described in the paper
for generating Fig 5.

# MaximalCapacity
It contains the code for reproducing S7.

This code can be modified to reproduce most of the figures on the paper.


