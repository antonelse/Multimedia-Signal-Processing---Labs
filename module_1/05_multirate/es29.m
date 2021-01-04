% Es 29

% Given the signal x(t)= A1 cos(2 pi f1 t) + A2 cos(2 pi f2 t):
% Create the signal x(n) as x(t) with t from 0 to 0.5 seconds, 
% sampled at Fs=8000 Hz. A1=0.7, A2=0.5, f1=1800 Hz, f2=3600 Hz
% Create the signal y(n) by resampling x(n) with 6000 Hz, without 
% using the MATLAB functions for automatic re-sampling. 
% Use N = 64 filter samples.
% Plot the magnitude of the DFTs of x(n), the upsampled signal, 
% the filtered signal and y(n) over 2048 samples 
% vs normalized frequency in [0, 1)

close all
clearvars
clc

%% parameters

A1 = .7;
A2 = .5;
f1 = 1800;
f2 = 3600;
Fs = 8000;
duration = 0.5;
time = 0:1/Fs:duration;

x = A1*cos(2*pi*f1*time) + A2*cos(2*pi*f2*time);

%% resampling without matlab functions

Fs_new = 6000;
[L, M] = rat(Fs_new / Fs); 

% first, upsampling
x_upsampled = zeros(1, length(x) * L);
x_upsampled(1:L:end) = x;

% filtering 
cutoff = min([1/(2*L), 1/(2*M)]);
cutoff_filter = 2*cutoff;
h = L*fir1(63, cutoff_filter);
x_f = filter(h, 1, x_upsampled);

% downsampling
y = x_f(1:M:end);

% %% resampling with MATLAB functions
% y_1 = decimate(interp(x, L), M);

%% DFTs

N = 2048;

Xf = fft(x, N);
Xup = fft(x_upsampled, N);
Xfilt = fft(x_f, N);
Yf = fft(y, N);

norm_freq_axis = 0:1/N:(N -1)/N;

figure; 
leg = {};
plot(norm_freq_axis, abs(Xf));
leg{1} = 'Original signal';
hold on;
plot(norm_freq_axis, abs(Xup));
leg{2} = 'Upsampled signal';
hold on;
plot(norm_freq_axis, abs(Xfilt));
leg{3} = 'Filtered signal';
hold on, 
plot(norm_freq_axis, abs(Yf));
leg{4} = 'Signal with Fs = 6000 Hz';
legend(leg);
grid
set(gca, 'fontsize', 16);



