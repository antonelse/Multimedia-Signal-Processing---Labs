% Es 24

% Write a MATLAB function typeOfFilter(b, a)
% that receives as input the numerator and denominator coefficients 
% of a causal filter H(z)=B(z)/A(z) and it returns: 
% -1 if the filter is not stable
% 1 if the filter is stable and it is minimum phase
% 0 if the filter is stable but it is not minimum phase
% If you test this function on a FIR filter, which is the output?
% Test the function on H(z) (see slides).

close all
clearvars
clc

%% define filter

B = [1, -2, 0, 2, -0.5, .2];
A = [1, 0.08, 0, 2];

%% allpass related filter

type = typeOfFilter(B, A);

% check
figure;
zplane(B, A);
