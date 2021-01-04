% Es 14

% Given H(z) = 1 - 0.8 z^-1
% Plot its impulse response

close all
clearvars
clc

%% define H(z)

% denominator coefficients (from a_0 to a_D)
A_z = 1;
% numerator coefficients (from b_0 to b_N)
B_z = [1, -.8];

%% H(f)

N = 100;
% you can use 'whole' to compute the Frequency response on 2pi 
[H_omega_2pi, omega_2pi] = freqz(B_z,A_z, N, 'whole');

figure,
plot(omega_2pi - pi, fftshift(abs(H_omega_2pi)));
xlabel('\omega');
title('amplitude');
grid

figure,
plot(omega_2pi - pi, fftshift(angle(H_omega_2pi)));
xlabel('\omega');
title('Phase');
grid

%% plot zeros and poles

% Be careful! If you first try computing the poles just giving 
% A_z to function 'roots', you will be wrong!
figure;
zplane(B_z, A_z);

%% function minmax phase

zeroes = roots(B_z);
is_minmax_phase(zeroes);







