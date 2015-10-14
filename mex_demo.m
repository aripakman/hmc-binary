% C++ implementation for Markov Random Fields (MRFs) of the 
% algorithm introduced in the NIPS 2013 paper "Auxiliary-variable exact 
% Hamiltonian Monte Carlo samplers for binary distributions" 
% by Ari Pakman and Liam Paninski. 
%

% The implementation uses the Gaussian augmentation.

% The log probability of the MRF is defined as 
%
% log p(S|r,M) = -sum_i r(i)*S(i) - sum_{i<j} M(i,j)*S(i)*S(j)  + const  
%
% with S(i) = +1 or -1. The sampler is called as 
% 
% [X,ll] = hmc_binary(M,r,L,P,last_y, seed);

% Input
% M:        d x d symmetric matrix of real numbers
% r:        d x 1 vector of real numbers
% L:        number of samples
% P:        the travel time of the particle is (P+.5)*pi
% last_y :  d x 1 vector of real numbers. See explanation below
% seed:     optional input for random seed. If absent, the seed is chosen from the C rand() function

% Output
% X:        d x L matrix, each column is a sample
% ll:       1 x L vector, with the log probability of each sample


%% last_y
% The Markov chain is defined over the continuos variables y, such that
% S = sign(y). The sampler requires initial values for y in the
% variable last_y. When the function returns, last_y will contain the last 
% values of y (last_y is passed by reference). 
% This is useful when the binary variables are part of a Gibbs sampling scheme.


%% Compile the source code for the mex function
mex  COMPFLAGS='$COMPFLAGS -Wall -std=c++11'  ...
     src_cpp/hmc_binary.cpp src_cpp/MRF.cpp ...
     src_cpp/HMC_BinarySampler.cpp src_cpp/BinaryDistribution.cpp
 

%% Example of use

d = 100;                    
beta = .7;
r = beta*(rand(d,1)-.5);
M = beta*(rand(d,d)-.5);
M = .5*(M+M');

L = 1e4;
P = 60;

last_y = ones(1,d);
[X,ll] = hmc_binary(M,r,L,P,last_y);

% plot samples
figure(9)
imagesc(X)


% plot magnetization and log probabilities for first 1000 samples

figure(10)
m = mean(X,1);
plot(m(1:1000))

figure(11)
plot(ll(1:1000))


