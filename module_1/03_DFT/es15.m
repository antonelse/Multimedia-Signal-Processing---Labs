% Es 15

% Given y(n) = -2 y(n-1) -  y(n-2) + x(n) + 2rho cos(theta)x(n-1) +
% rho^2 x(n-2)
% rho = 0.9, theta = pi/8
% The sequence is defined for n in [0, 1e4].
% Which is the expression of h(n)? 
% Which is the amplitude of H(f)? 
% Compute the DFT of h(n) using the matrix product.

% which kind of filter is it? (Low pass, band pass, high pass...)
% Is it a FIR or IIR filter?

close all
clearvars
clc

%% parameters

rho = 0.9;
theta = pi/8;

%% define H(z)

% denominator coefficients (from a_0 to a_D)
A_z = [1, 2, 1];
% what would happen if A_z = 1? Try by yourself.
% A_z = 1;
% numerator coefficients (from b_0 to b_N)
B_z = [1, 2*rho*cos(theta), rho^2];

%% h(n)

N = 1e4;
n = 0:N-1;
delta = zeros(size(n));
delta(1) = 1;
h = filter(B_z, A_z, delta);

%% amplitude of H(f)

% you can use 'whole' to compute the Frequency response on 2pi 
[H_f_2pi, omega_2pi] = freqz(B_z, A_z, N, 'whole');

figure, 
% plot(omega_2pi, (abs(H_f_2pi)));
% in order to better see the behaviour,you can use 'semilogy'
semilogy(omega_2pi, (abs(H_f_2pi)))
xlabel('\omega');
ylabel('|H(\omega)|');
set(gca, 'fontsize', 16);
grid

%% compute the DFT of h using matrix product

% just with one line 
W = exp(-1i * 2 * pi /N.*n'.*n);
H_k = W * h';

%% plot between 0 and 2pi.

hold on;
semilogy(2*pi*n/N, abs(H_k), '--')

%% plot between -pi, pi

figure; 
semilogy(omega_2pi - pi, abs(fftshift(H_f_2pi)))
hold on; 
semilogy(2*pi*n/N - pi, fftshift(abs((H_k))), '--');

%% compute H_k using fft

H_fft = fft(h);
hold on;
semilogy(2*pi*n/N - pi, fftshift(abs((H_fft))),  '--o');



