clear;clc;close all;
load('true_param.mat')
load("ref_traj.mat")

P = traj(1,:,1);
t = zeros(2,60000);

for i=1:60000
    t(:,i) = traj(1,:,i)' - A(:,:,i)*P';
end
