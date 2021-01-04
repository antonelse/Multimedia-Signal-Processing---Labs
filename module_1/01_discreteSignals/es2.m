% Es 2 from 7/11/2018

% Given the signal x(t) = Acos(2pift)
% 1. Write the script es2.m to create the signal x(n) as x(t)
% from 0 to 0.5 seconds, sampled at Fs (sampling rate) =
% 1000Hz; A = 0.8, f = 50Hz.
% 2. Write the function sinusoid.m which takes as input the
% time-axis, the amplitude, the frequency, the phase of a
% discrete sinusoid and return the signal.
% 3. Generate the same signal as 1. with sinusoid.m
% 4. In es2.m, plot the signal as a function of n.
% 5. In es2.m, plot the signal as a function of time samples.

close all
clearvars
clc

%% define parameters

A = 0.8;
f = 50;
Fs = 1000;
tmax = 0.5;
tmin = 0;
time = tmin:1/Fs:tmax;

x = A * cos(2*pi*f*time);

%% generate the signal with function 'sinusoid.m'

x1 = sinusoid(time, A, f, 0);

%% plot

figure(1);
plot(x);
hold on; 
plot(x1, '--')
xlabel('n (sample index)');
title('Sinusoid');

figure(2);
plot(time, x);
hold on; 
plot(time, x1, '--')
xlabel('time [seconds]');
title('Sinusoid');





