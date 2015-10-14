hmc-binary
=======

MATLAB and C++ implementations of the exact HMC sampler for binary distributions using a Gaussian augmentation. The algorithm was introduced in the NIPS 2013 paper "Auxiliary-variable exact Hamiltonian Monte Carlo samplers for binary distributions" by Ari Pakman and Liam Paninski. 

Paper available at http://arxiv.org/abs/1311.2166

The script ising_1d_demo.m runs the MATLAB implementation and reproduces the 1D Ising model example in the paper.

The script mex_demo.m compiles and runs the C++ version for binary Markov Random Fields (Boltzmann machines). For other binary distributions, one can write a class derived from HMC_BinarySampler (similar to the MRF class). 



