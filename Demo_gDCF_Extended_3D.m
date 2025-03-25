close all; clear all; clc

%Load in sample trajectory
%[Nsample,Nshot,Ndim,Nt] =size(k) %Ndim = 3 for x,y,z,Nt= timeframe
% k space trajectory normalized between [-0.5,0.5,0.5]

data=load('.\Cones_106Shots.mat'); % enter filepath
k_space_traj = data.proj_grid; 
size(k_space_traj)
%recon matrix size
N = 128;

%compute DCF
gDCF = gDCF_extended_3D(k_space_traj,N);

%Plot DCF for one time frame
figure; plot3( k_space_traj(:,:,1), k_space_traj(:,:,2), k_space_traj(:,:,3)) 
title(['k-space Traj'])

figure(2);
plot(gDCF(:,1,1),'r','LineWidth',2)
title(['gDCF-ext'])