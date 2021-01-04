% Es 12

% Given y(n) = x(n) - a x(n-1) + b y(n-1) 
% Which is the expression of H(z)?
% Which is the value of h(0)?
% Compute and plot the zeros and poles.
% Plot h(n) for n in [0, 50] for a = 0.5 and b = 0.2. 
% Plot h(n) for n in [0, 50] for a = 1.2 and b = 0.2
% Plot h(n) for n in [0, 50] for a = 1.2 and b = 1.1
% In which situations is the system stable?

close all
clearvars
clc

%% parameters

a = 1.2;
b = 1.1;

% denominator
A_z = [1, -b];
% numerator
B_z = [1, -a];

%%% h(0) must be 1.

%% zeros and poles

% be careful in defining the variable for the zeros with a name different 
% from 'zeros', which is a MATLAB function.
zeroes = roots(B_z);
poles = roots(A_z);

%% plot them

figure; 
zplane(zeroes, poles);
grid

%% h(n) 

n = 0:50;
delta = zeros(size(n));
delta(1) = 1;
h = filter(B_z, A_z, delta);

%% plot h(n)

figure, 
stem(n, h);






