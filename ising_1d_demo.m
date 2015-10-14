% Author: Ari Pakman
% 
% Demo for the augmented-variable HMC sampler introduced in the NIPS 2013 paper 
% "Auxiliary-variable exact Hamiltonian Monte Carlo samplers for binary distributions" by Ari Pakman and Liam Paninski 
%
% This MATLAB implementation reproduces one of the examples of the paper.
% See also the C++ implementation with mex function.

clear
addpath('src_matlab')

d =400;
temp =.45;  % temperature for the 1D Ising model

is1 = Ising1D(d,temp);   % create 1D Ising object 

t =12.5;        % computational cost of every sample. useful to compare HMC vs Metropolis
L=500;          % number of samples 


% run the Metropolis sampler
tic
[IsGs, IsGL1]  = MetroGibbs_binary(is1,L, d*t);
toc

T=t*pi;  % time to run each HMC iteration

% run the HMC sampler (with Gaussian augmentation)
tic
[IsHs, IsHL1] = HMC_binary(is1,t*pi,L);
toc


%% Plot the results 

maG = mean(IsGs,1);
maH = mean(IsHs,1);

fig=figure(82);
clf

bwmap = [.97 .97 1 ; 0 0 0];
colormap(bwmap);

subplot(411)
hold on
plot(maG(1:end))
plot(maH(1:end),'r')
title('Magnetization');
grid
box on

subplot(412)
cla
hold on
plot(-IsGL1(1:end))
plot(-IsHL1(1:end),'r')
title('Energy');
leg=legend('Metropolis', 'HMC');
set(leg, 'FontSize', 7)
grid
box on

subplot(413)
imagesc(IsGs(:,1:end))
title('Metropolis');

subplot(414)
imagesc(IsHs)
ti=title('HMC');
xlabel('Iteration')

