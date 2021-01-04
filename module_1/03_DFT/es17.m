% Es 17

% Given a sinusoidal signal with frequency = 2Hz and duration = 1.3 seconds, sampled with Fs = 50 Hz
% Plot the first N = 50 samples of the sinosoid vs time.
% Compute its FFT on N samples and visualize it vs frequencies.
% Compute the FFT of the complete signal and visualize it.
% Do the two Fourier spectra coincide? 
% Which frequencies correspond to peaks in Fourier?  Why?
% How to find again the same Fourier peak in the second spectrum? 

close all
clearvars
clc

%% parameters

f0 = 2;
Fs = 50;
duration = 1.3; 
time = 0:1/Fs:duration;

%% define signals

s = cos(2*pi*f0*time);

%% plot the signal on the first 50 samples

N = 50;

figure;
plot(time(1:N), s(1:N));

%% FFT on the first 50 samples.

S_f = fft(s(1:N));
freq = 0:1/N*Fs: Fs - 1/N*Fs;

figure;
stem(freq - Fs/2, fftshift(abs(S_f)));

%% FFT of the complete signal

S_complete = fft(s);
N_complete = length(s);
freq_complete = 0:1/N_complete * Fs: Fs - 1/N_complete * Fs;

figure;
stem(freq_complete - Fs/2, fftshift(abs(S_complete)));

%% let's plot the complete signal

figure; 
plot(time, s);

%% how to find the actual peak in Fourier?

% zero padding helps finding the actual frequency of the sinusoid.
% since you are interpolating samples... You can find the right peak
% the final signal is just the convolution of the peaks due to the cosine
% by the periodic sinc.

N_tot = 1000;
num_pad = N_tot - N_complete;
s_pad = padarray(s, [0, num_pad], 'post');
S_pad = fft(s_pad);

N_pad = N_complete + num_pad;
freq_pad = 0:1/N_pad * Fs: Fs - 1/N_pad * Fs;

figure; 
stem(freq_pad - Fs/2, fftshift(abs(S_pad)));
