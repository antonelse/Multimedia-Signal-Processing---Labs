% Es 16

% Given x(n) = rectangular pulse of width L, n in [0, L-1]
% Given y(n) = x(n)
% Compute z(n) as the linear convolution between x and y.
% Compute the cyclic convolution between x and y

close all
clearvars
clc

%% define signals

L = 6;
x = ones(1, L);
n_x = 0:L-1;
y = x;
n_y = 0:L-1;

%% convolve them linearly

z = conv(x, y);
n_z = 0:n_x(end) + n_y(end);

figure, stem(n_z, z)

%% cyclic convolution

z_c = cconv(x, y, L);

figure;
stem(0:L-1, z_c, '--')

%% fft

X = fft(x);
Y = fft(y);

Z = X.*Y;

z = ifft(Z);

hold on,
stem(0:L-1, z)

%% compute the cyclic convolution as a linear convolution.

N = L + L -1;

% padding the array with zeros: do not add anything on rows and add (N-L)
% columns on the right.
x_pad = padarray(x, [0, N-L], 'post');
y_pad = padarray(x, [0, N-L], 'post'); 
z_pad = cconv(x_pad, y_pad);
n_pad = 0:N-1;

figure, stem(n_z, z_pad(n_pad >=0 & n_pad <N));

%% with DFT

X_pad = fft(x_pad);
Y_pad = fft(y_pad);
Z_pad = X_pad.*Y_pad;

z_pad_ifft = ifft(Z_pad);

hold on, stem(n_z, z_pad_ifft(n_pad >=0 & n_pad <N), '--');


