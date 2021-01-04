% Es 23

% Write a MATLAB function allpass.m which has the form:
% [z_out, p_out, b_out, a_out] = allpass(b,a)
% Input: b, a = numerator and denominator of H(z)
% Output: z_out, p_out, b_out, a_out = zeros, poles, numerator, 
% denominator of the allpass transfer function related to H(z)
% Use the function allpass.m to compute the allpass transfer function related to the causal filter
% Plot the magnitude response of the filter vs normalized frequencies using N = 512 samples
% How do you expect the phase to behave?

close all
clearvars
clc

%% define filter

B1 = [1, 3];
A1 = [1, .5];

%% allpass related filter

[z_ap, p_ap, B_ap, A_ap] = allpass(B1, A1);

%% plot the magnitude response as a function of omega

N = 512;
[Hap, omega] = freqz(B_ap, A_ap, N, 'whole');

figure;
plot(omega./(2*pi), abs(Hap));

% phase behaviour: a causal stable allpass system has always maximum phase
% zeros (they are the reciprocal complex conjugate of the poles)
% therefore the phase has jumps = 2*pi
figure;
plot(omega./(2*pi), angle(Hap));

% %% plot the zplane
% 
% figure;
% zplane(B_ap, A_ap)
