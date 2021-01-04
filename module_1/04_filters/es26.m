% Es 26

% Given x as a cosine wave sampled at Fs = 8KHz, duration 1 second, 
% amplitude 1.5, frequency 1.1KHz, phase 45 deg.
% Plot the sinusoid vs time
% Compute y as x filtered with a low-pass filter 
% with normalized cut-off frequency of 0.4 and 64 weights
% Apply a Hanning window to select the first 512 samples of y
% Plot the magnitude of the DFT of the windowed y versus frequency in Hz.
% If you change the cut-off frequency to 0.05, 
% what do you expect to see in the spectrum of y?

close all
clearvars
clc

%% parameters

A = 1.25;
% convert the phase in radians
phase = 45*pi/180;
f0 = 1.1e3;
Fs = 8e3;
duration = 1;
time = 0:1/Fs:duration;

%% signal

x = A*cos(2*pi*f0*time + phase);

% % fft magnitude of the signal
% % you can use plot or semilogy
% figure, semilogy([0:1/length(x) *Fs:(length(x)-1)/length(x)*Fs] - Fs/2, abs(fftshift(fft(x))));

%% Compute y as x filtered with a low-pass filter with normalized cut-off frequecy of 0.4 and 64 weights 4

N_filter = 64;
cutoff = .4;
% NB: to obtain the cutoff of the filter, multiply by 2
cutoff_filter = cutoff * 2;

% create the FIR filter
h = fir1(N_filter -1, cutoff_filter);

% % frequency response of the filter vs normalized frequencies
% [H, omega] = freqz(h, 1, N_filter, 'whole');
% figure, semilogy((omega-pi)./(2*pi), fftshift(abs(H)));

% filtered signal
y = filter(h, 1, x);

% % frequency response of the output signal vs frequency [Hz]
% figure, semilogy([0:1/length(y)*Fs:(length(y)-1)/length(y)*Fs] - Fs/2, abs(fftshift(fft(y))));

%% hanning window to select the first 512 samples

Nw = 512;
w = window(@hanning, Nw);
% try also rectwin to see what changes
y_w = y(1:Nw) .* w';

%% Plot the magnitude of the DFT of the windowed y vs Hz

Yw = fft(y_w);
figure, plot([0:Fs/Nw:(Nw-1)/Nw *Fs] - Fs/2, abs(fftshift(Yw)));




