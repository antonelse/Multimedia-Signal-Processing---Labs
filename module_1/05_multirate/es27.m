% Es 27

% Given x(n) defined as the sum of two sinusoidal signals, 
% sampled at Fs = 500 Hz with duration 3 seconds, one with frequency 50 Hz 
% and the other one with frequency 100Hz:
% Downsample x(n) with downsampling factor M = 4
% Decimate x(n) with decimation factor M = 4, using a FIR filter with order 64. 
% Plot the DFTs of x(n), of the downsampled and of the decimated signals 
% vs frequency [Hz] in the same figure and comment on the results.
% Try also M = 2 and see what happens

close all
clearvars
clc

%% parameters

Fs = 500;
f_0 = 50;
f_1 = 100;
duration = 3;
time = 0:1/Fs:3; % - 1/Fs; subtract 1/Fs to see a perfect delta in the DFT
x = cos(2*pi*f_0 * time) + cos(2*pi*f_1 * time);

%% downsample x with a factor M

M = 4;
x_downsampled = x(1:M:end);

%% decimate the signal with a factor M 

% consider a filtering with cutoff = Fs/(2*M) 
% --> normalized cut-off frequency = 1/(2*M)
% --> cut-off frequency for the filter = 1/M

lpf = fir1(64, 1/M);
% filtered signal
x_f = filter(lpf, 1, x);
% downsample the signal
x_decimated = x_f(1:M:end);

%% DFT of x

Xf = fft(x);
N = length(x);
freq_axis = 0:Fs/N:Fs*(N-1)/N;

figure;
plot(freq_axis, abs(Xf));
leg{1} = 'Original signal';

%% DFT of downsampled and decimated signals

N_new = length(x_downsampled);

Xf_downsampled = fft(x_downsampled);
Xf_decimated = fft(x_decimated);

freq_axis = 0:Fs/N_new : Fs*(N_new -1)/N_new;

hold on;
plot(freq_axis, abs(Xf_downsampled), '--x');
leg{2} = 'Downsampled signal';
hold on;
plot(freq_axis, abs(Xf_decimated), '--o');
leg{3} = 'Decimated signal';
legend(leg);








