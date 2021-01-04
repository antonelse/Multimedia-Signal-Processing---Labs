% Es 16

% Given a FIR filter h(n) = 1, n in [0, 19].
% Evaluate the DFT of the filter
% Visualize the DFT versus normalized frequency [0.5, 0.5).
% Pad the array with zeros until reaching 100 samples.
% Evaluate the DFT of the padded h(n) and visualize it.
% Are the two DFTs equal? Comment on the results.

close all
clearvars
clc

%% define signals

h = ones(1, 20);
N = length(h);
n_h = 0:N-1;

%% fft result

H_f = fft(h);

figure, stem((n_h)./N - .5, fftshift(abs(H_f)))

%% padding with zeros

num_zeros = 1000 - N;
h_pad = padarray(h, [0, num_zeros], 'post');
N_pad = length(h_pad);

H_pad = fft((h_pad));

figure, stem([0:(N_pad - 1)]./(N_pad) - .5, fftshift(abs(H_pad)));

% padding with zeros = convolve in frequency domain by a periodic sinc...
% We are interpolating samples, increasing the density of Fourier spectrum

