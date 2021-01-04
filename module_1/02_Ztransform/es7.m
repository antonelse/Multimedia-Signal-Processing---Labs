% Es 7

% Generate g(n) as a set of 10000 realizations of random variables distrubted as ~ N(0, 1).
% Compute the mean and the variance
% Create h(n) = a * g(n) + b, with a = 0.5; b = 4. 
% Is h(n) still a Gaussian random variable?
% Compute the mean and the variance of h(n).
close all
clearvars
clc

%% parameters

N = 1e4;
a = .1;
b = 4;

%% create g(n)

g = randn(1, N);

%% estimate mu and var

mu_g_est = mean(g);
var_g_est = var(g);

%% plot an estimate of the pdf

figure; 
h_g = histogram(g, 'Normalization', 'pdf');
% verify that pdf sums up to 1
pdf_area = sum(h_g.Values)*h_g.BinWidth; 
title(sprintf('Gaussian pdf, estimated mu = %f and std = %f', mu_g_est, sqrt(var_g_est)));
set(gca, 'fontsize', 14)

%% create h(n)

h = a * g + b;
% h is gaussian as well! If operations are linear, everything remains gaussian

% theoretical mu and var: mu = 4, var = .01;
mu_h_est = mean(h);
var_h_est = var(h);
std_h_est = std(h); 

%% plot an estimate of the pdf of h

figure; 
h_h = histogram(h, 'Normalization', 'pdf');
title(sprintf('PDF of h, estimated mu = %f and std = %f', mu_h_est, sqrt(var_h_est)));
set(gca, 'fontsize', 14)

