% Es 5

% Generate the signal x(n) = u(n-5) - u(n-10), considering n = 1:15.
% Generate the periodic signal xp(n) with period N = 15, considering n = 1:200.
% Hint: Consider using repmat instead of for loops.
% Plot the periodic signal xp(n) considering only 8 periods.

close all
clearvars
clc

%% define one period of the signal

N = 15;
n = 1:15;
x = zeros(1, N);
x(n>= 5 & n <10) = 1;

%% plot the signal with the function stem.

figure, 
stem(x)
grid

%% create the periodic signal

n_max = 200;
% how many signal periods until n_max?
N_p = floor(n_max / N); 

% create the periodic signal until N_p + 1 periods (because 200 exceeds N_p*N)
% repmat(x, #repetitions on rows, #repetitions on columns)
x_p = repmat(x, 1, N_p + 1);

% check the length of x_p: it exceeds 200
l_p = length(x_p);

% cut the signal at n = n_max
x_p = x_p(1:n_max);

%% plot the periodic signal xp(n) considering only 8 periods

figure, 
stem(1:N*8, x_p(1:N*8))
grid







