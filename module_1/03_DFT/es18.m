% Es 18
% Given a sinusoidal signal with frequency = 2Hz and duration = 1 seconds, 
% sampled with Fs = 50 Hz
% Add this signal to a second sinosoid with frequency 2.2Hz 
% and the same durantion and sampling rate.
% Plot the global signal vs time.
% Compute the FFT. Can you see the two peaks?
% How to see the correct Fourier spectrum?

close all
clearvars
clc

%% parameters

f0 = 2;
f1 = 2.2;
Fs = 50;
N = 50; 
time = 0:1/Fs:1/Fs*(N-1);

%% define signals

s0 = cos(2*pi*f0*time);
s1 = cos(2*pi*f1*time);
s = s0 + s1;

%% plot the signal

figure;
plot(time, s);

%% FFT of s

S_f = fft(s);
freq = 0:1/N*Fs:Fs - 1/N*Fs;

figure;
stem(freq - Fs/2, fftshift(abs(S_f)));

%% zero padding?

% we can increase the density in Fourier domain, but we will not see the
% two peaks because the resolution is fixed by the time window!
% zero-padding does not modify the size of the rectangular window.

N_tot = 1000;
num_pad = N_tot - N;
s_pad = padarray(s(1:50), [0, num_pad], 'post');
S_pad = fft(s_pad);

freq_pad = 0:1/N_tot * Fs: Fs - 1/N_tot * Fs;

figure; 
stem(freq_pad - Fs/2, fftshift(abs(S_pad)));

%% increase the number of measurements

% increasing measurements modifies the size of the rectangular
% window --> frequency resolution enhances.
% to see the 2 exact peaks: 
% you should have enough resolution
% the number of samples must be a multiple of the two periods

N = 250;
time = 0:1/Fs:1/Fs*(N-1);

s0 = cos(2*pi*f0*time);
s1 = cos(2*pi*f1*time);
s = s0 + s1;

S_f = fft(s);
freq = 0:1/N*Fs:Fs - 1/N*Fs;

figure;
stem(freq - Fs/2, fftshift(abs(S_f)));





