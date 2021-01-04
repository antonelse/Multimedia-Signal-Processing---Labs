% Es 28

% Given the downsampled signal defined in Es27 with M = 4, 
% Create the signal x1 by upsampling the signal with a factor L = 4
% Given the decimated signal defined in Es27 with M = 4,
% Create the signal x2 by interpolating the signal with a factor L = 4, 
% using a FIR filter with order 64.
% Open a figure and create three subplots:
% In 1° subplot, plot the stem of the original signal x(n) until N = 120 
% time samples, x-axis in seconds.
% In 2° subplot, plot the stem of the downsampled and decimated signals 
% with the same temporal duration as above
% In 3° subplot, plot the stem of x1 and x2 with the same temporal duration


close all
clearvars
clc

%% parameters

Fs = 500;
f_0 = 50;
f_1 = 100;
duration = 3;
time = 0:1/Fs:3;
x = cos(2*pi*f_0 * time) + cos(2*pi*f_1 * time);

%% downsample x with a factor M

M = 4;
x_downsampled = x(1:M:end);

%% upsample the signal x_downsampled by a factor L = 4

L = 4;
x1 = zeros(1, length(x_downsampled) * L);
x1(1:L:end) = x_downsampled;

%% decimate the signal with a factor M 

% consider a filtering with cutoff = Fs/(2*M) 
% --> normalized cut-off frequency = 1/(2*M)
% --> cut-off frequency for the filter = 1/M

lpf = fir1(64, 1/M);
% filtered signal
x_f = filter(lpf, 1, x);
% downsample the signal
x_decimated = x_f(1:M:end);

%% interpolate the signal x_decimated by a factor L

x2 = zeros(1, length(x_decimated) * L);
x2(1:L:end) = x_decimated;

% consider a filtering with cutoff = 1/(2*L) 
% --> cut-off frequency for the filter = 1/L

lpf = fir1(64, 1/L);
% remember to put the gain L in interpolation
x2 = L*filter(lpf, 1, x2);

%% plot the signals

N = 120;

figure;
leg = {};
% three rows, 1 column, subplot idx
subplot(3, 1, 1);
stem(time(1:N), x(1:N));
leg{1} = 'Original signal';
legend(leg);

subplot(3, 1, 2);
leg = {};
% for plotting the downsampled signal using the same time axis as the
% original one, define the equivalent sampling frequency
Fs_down = Fs/M;
% time axis goes from 0 until time(N)
time_down = 0:1/Fs_down:time(N);
stem(time_down, x_downsampled(1:length(time_down)));
leg{1} = 'Downsampled signal';
hold on;
% decimated signal
stem(time_down, x_decimated(1:length(time_down)), '--*');
leg{2} = 'Decimated signal';
legend(leg);

subplot(3, 1, 3)
leg = {};
% upsampled signal x1
stem(time(1:N), x1(1:N));
leg{1} = 'Upsampling the downsampled signal';
hold on,
% interpolated signal x2
stem(time(1:N), x2(1:N), '--*');
leg{2} = 'Interpolating the decimated signal';

% by interpolating the decimated signal, we can exactly reconstruct 
% the sinusoid with frequency < min(pi/M, pi/L)
% NB: in the reconstruction, there is a delay introduced by the filter
% which cannot be centered in 0, as the filter must be causal. 
% the total delay is given by the sum of the two delays introduced by the
% FIR filter used in decimation and the FIR filter used in interpolation.
% if the filter has 65 taps --> the filter peak is centered in sample #33
% total delay = 33 samples for decimation + 33 samples for interpolation 
% --> 66 samples
% to find the delay with MATLAB, use the function max.
% the second output finds the index corresponding to the maximum value of
% the array.
[~, filter_delay] = max(lpf);
total_delay = 2*filter_delay;
hold on; 
stem(time(1:N), cos(2*pi*f_0 * time(total_delay + 1:total_delay +N)));
leg{3} = 'Band-limited signal, f = f_0, delayed by 66 samples';
legend(leg);

%% DFT of the signals

figure;
leg = {};
freq_axis = 0:1/N:(N-1)/N;
% three rows, 1 column, subplot idx
subplot(3, 1, 1);
plot(freq_axis, abs(fft(x(1:N))));
leg{1} = 'Original signal';
legend(leg);

subplot(3, 1, 2);
leg = {};
% for plotting the downsampled signal using the same time axis as the
% original one, define the equivalent sampling frequency
Fs_down = Fs/M;
% time axis goes from 0 until time(N)
time_down = 0:1/Fs_down:time(N);
freq_axis_down = 0:1/length(time_down):(length(time_down)-1)/length(time_down);
plot(freq_axis_down, abs(fft(x_downsampled(1:length(time_down)))));
leg{1} = 'Downsampled signal';
hold on;
% decimated signal
plot(freq_axis_down, abs(fft(x_decimated(1:length(time_down)))), '--*');
leg{2} = 'Decimated signal';
legend(leg);

subplot(3, 1, 3)
leg = {};
% upsampled signal x1
plot(freq_axis, abs(fft(x1(1:N))));
leg{1} = 'Upsampling the downsampled signal';
hold on,
% interpolated signal x2
plot(freq_axis, abs(fft(x2(1:N))), '--*');
leg{2} = 'Interpolating the decimated signal';
legend(leg);






