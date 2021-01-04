% Es 9

% Given x(n) = [3, 11, 7, 0, -1, 4, 2] , n in [-3, 3]
% Create y(n) = x(n - 5), n in [0, 10], without using circshift or for loops.
% Create y(n)= 1/3 sum[m = 0, 1, 2] x(n-m)
% Hint: y(n) has the form of a convolution...

close all
clearvars
clc

%% signals

x = [3, 11, 7, 0, -1, 4, 2]; 
n_x = [-3:3];

%% y(n) = x(n-5)

n_h = 0:10;
delta_5 = zeros(size(n_h));
delta_5(6) = 1;

% support of the convolution:
n_conv = n_x(1) + n_h(1) : n_x(end) + n_h(end);

y = conv(delta_5, x);
y = y(n_conv>=0 & n_conv <= 10);

stem(0:10, y);

%% y(n) = sum...

h = zeros(size(n_h));
h(1:3) = 1/3;

y = conv(x, h);
y = y(n_conv>=0 & n_conv <= 10);

figure; 
stem(n_x, x);
figure; 
stem(0:10, y);

