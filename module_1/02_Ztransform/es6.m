% Es 6

% Generate noise(n) as a set of uniformly distributed random variables 
% between -A and A, A=0.05. The time axis has a 2 seconds duration with Fs=11.025 Hz.
% Generate x(n) as sinusoidal sequence with frequency 220 Hz
% Generate the signal y(n) = x(n) + noise(n) 
% Normalize the signal in [-1, 1]
% Play y(n), testing different values of A
close all
clearvars
clc

%% parameters

Fs = 11025;
A = 10; 0.05;
duration = 2;
time = 0:1/Fs:2;
f0 = 220;

%% generate the noise

noise_01 = rand(size(time));
% how to change the dynamics from [0, 1] to [-A, A]?
% to pass from [0, 1] to [a, b]-->
% MULTIPLY BY THE TOTAL RANGE (b-a) AND SUM a
noise = noise_01*2*A - A;

%% check the histogram.

figure;
histogram(noise);

%% generate the signal x

x = cos(2*pi*f0*time);

%% generate the signal y

y = x + noise;

%% normalize the signal in [-1, 1]

% how to change the dynamics?
% to pass from [min(y), max(y)] to [-1, 1]:
% it's always better to convert dynamics to [0, 1] -->
% SUBTRACT min(y)and DIVIDE BY (max(y) - min(y)) 
% then convert to [-1, 1] -->
% MULTIPLY BY THE TOTAL RANGE (2) AND SUBTRACT -1

% GENERAL RULE FOR SCALING THE DYNAMICS
% y_out = (y_in - min(y_in)) / (max(y_in) - min(y_in)) * 
% (max(y_out) - min(y_out)) + min(y_out)

y_n = ((y - min(y))/(max(y)- min(y))) * 2 - 1;

%% plot

figure; 
plot(time(1:200), y_n(1:200));
ylabel('sinusoid + uniformly distributed random noise');
xlabel('time');
title(sprintf('A = %f', A));

%% play the sinusoids

sound(y_n,Fs);
