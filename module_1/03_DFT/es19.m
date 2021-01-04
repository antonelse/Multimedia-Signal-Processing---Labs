% Es 17

% circular convolution.

% Given x(n) = delta(n-2), defined for n in [0, 4]
% Given y(n) = [5, 4, 3, 2, 1], n in [0, 4].
% Compute z(n) as the linear convolution between x and y.
% Which is the support of z(n)?
% Compute the cyclic convolution between x and y.
% Try using matlab function cconv and check the result.


close all
clearvars
clc

%% define signals

x = [0, 0, 1, 0, 0];
n_x = 0:4;
y = [5, 4, 3, 2, 1];
n_y = 0:4;
N = length(x);

%% convolve them linearly

z = conv(x, y);
n_z = n_x(1) + n_y(1):n_x(end) + n_y(end);

figure, stem(n_z, z)

%% cyclic convolution

% let us create the periodic extension of y (one period on the left and one
% one the right it's enough)

% repmat: dimension of y_tilde should be: 1*num_rows(y), 3*num_cols(y)
y_tilde = repmat(y, 1, 3);
% support of the periodic extension
n_tilde = -N:2*N-1;

% now compute linear convolution
z_p = conv(x, y_tilde);
n_zp = n_tilde(1)+n_x(1):n_tilde(end) + n_x(end);

% take the cyclic convolution between 0 and N-1
z_c = z_p(n_zp >= 0 & n_zp <N);

figure, stem(0:N-1, z_c);

%% matlab function

z_mat = cconv(x, y, N);
hold on, 
stem(0:N-1, z_mat, '--')







