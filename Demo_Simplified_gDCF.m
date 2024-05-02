close all; clear all; clc

%Load in sample spiral trajectory
%data size = [no of readout, no of spiral arms, time frame]
%trajectory is normalized within [-0.5 0.5]

data=load('...\k-space traj.mat'); % enter filepath
k_space_traj = data.traj; 

%recon matrix size
N = 128;

%compute DCF
gDCF = Simplified_gDCF(k_space_traj,N);

%Plot DCF for one spiral arm in timeframe #1
figure;
plot(gDCF(:,1,1), 'r','LineWidth',2,'DisplayName', 'gDCF-ext')