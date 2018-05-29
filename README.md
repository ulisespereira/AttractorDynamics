# Attractor dynamics for networks with learning rules inferred from in vivo data (U. Pereira & N. Brunel, Neuron 2018)

This repository contains three directories:

1. Dynamics
2. CapacityFung
3. MaximalCapacity
 
With a few modifications, this code reproduces most of the figures on the paper. Send me an 
e-mail if you have further inquiries.

# Dynamics
It contains the code for reproducing Fig 6 A-F. 
The results on Fig 3 and Fig 4 A are reproduced by changing the value of the parameter A from 9.5 to 3.5 and setting the input current during the delay period as:

1. Standard normal iid for a novel stimulus
2. One of the stored patterns for familiar stimulus.


For generating fig4.pdf, first run `python main.py` on the `numeric/network_simulation` directory. The simulation
will take a couple of hours depending on your machine. It will create two large files `the_dynamics.p' (~800MB) and `the_overlaps.p' (~1.15GB). They store the dynamics of 100 neurons and all the overlaps for 10 realizations of the 8s dynamics.


# Capacity_vs_g
It contains the code for reproducing Fig 4 B-C. 
The same code can be used to perform the parameter exploration in Fig 5.

# MaximalCapacity
It contains the code for reproducing S7.



