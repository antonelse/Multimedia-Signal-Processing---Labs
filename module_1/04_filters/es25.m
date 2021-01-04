% Es 25

% Given the filter with 
% B(z)=[1, -1.98, 1.77, -0.17, 0.21, 0.34], 
% A(z)=[1, 0.08, 0.40 ,0.27]
% Compute the allpass-minimum phase decomposition of H(z)
% Check the results using ?zplane?
% Plot the magnitude of Hap(omega) using N = 1024 samples
% Plot the first 50 samples of h_min(n)

close all
clearvars
clc

%% define filter

B = [1, -1.98, 1.77, -0.17, 0.21, 0.34];
A = [1, 0.08, 0.40 ,0.27];

%% all-pass minimum phase decomposition

% compute the zeroes and the poles
zeroes = roots(B);
poles = roots(A);

% The minimum phase system contains zeroes and poles inside the unit circle
% plus the conj reciprocal of zeros outside the unit circle
z_min = zeroes(abs(zeroes) <= 1);
p_min = poles(abs(poles) <=1);

% check if there are zeros outside the unit circle
if any(abs(zeroes) > 1)
    z_min = [z_min; 1./conj(zeroes(abs(zeroes) > 1))];
    % allpass decomposition
    p_ap = 1./conj(zeroes(abs(zeroes) > 1));
    z_ap = zeroes(abs(zeroes) > 1);
end

B_min = poly(z_min);
A_min = poly(p_min);

B_ap = poly(z_ap);
A_ap = poly(p_ap);

%% check on zplane

figure;
zplane(B_min, A_min);

figure;
zplane(B_ap, A_ap);

%% magnitude of Hap

N = 1024;
[Hap, omega] = freqz(B_ap, A_ap, N, 'whole');
figure, plot(omega, abs(Hap));

%% h min

% minimum phase systems have always the energy concentrated in the first
% temporal samples
N = 50;
delta = zeros(1, N);
delta(1) = 1;
hmin = filter(B_min, A_min, delta);
figure, stem(hmin)








