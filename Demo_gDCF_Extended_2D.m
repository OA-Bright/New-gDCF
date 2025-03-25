close all; clear all; clc

%Load in sample trajectory
%data size = [Nreadout, Nshots, Ntime_frame]
%trajectory is normalized within [-0.5 0.5]

data=load('.\Radial_2D_36_spokes_GA.mat'); % enter filepath
k_space_traj = data.k_rad; 
size(k_space_traj)
%recon matrix size
N = 384;

%compute DCF
gDCF = gDCF_extended_2D(k_space_traj,N);

%Plot Traj and DCF for one timeframe 
figure(1)
plot(k_space_traj(:,:,1))
title(['k-space Traj'])

figure(2);
plot(gDCF(:,1,1),'r','LineWidth',2)
title(['gDCF-ext'])