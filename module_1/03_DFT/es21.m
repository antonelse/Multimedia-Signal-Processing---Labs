% Es 21

% overlap and add and overlap and save methods. 
% Given x(n) = n + 1, n in [0, 19]
% h(n) is a FIR filter, h(n) = [1, 0, -1], n in [0, 2]
% Compute the linear convolution between x and h
% Compute the same result with the overlap and add method, using L = 6.
% Compute the same result with the overlap and save method, using L = 6.

close all
clearvars
clc

%% define signals

n_x = 0:20;
x = n_x + 1;
h = [1, 0, -1];
n_h = 0:2;

%% linear convolution

y = conv(x, h);
n_y = n_x(1) + n_h(1) : n_x(end) + n_h(end);

figure, stem(n_y, y);

%% overlap and add method.

L = 6;
Lconv = L + length(h) - 1;

% how many non-overlapping blocks of length 6 are there in x? 
num_blocks = ceil(length(x) / L);

if length(x) < num_blocks * L
    % pad x with zeros 
    x_pad = padarray(x, [0, num_blocks*L - length(x)], 'post');
else
    x_pad = x;
end

y_oa = zeros(1, Lconv*num_blocks);
% you can directly tell matlab the fft samples. In case the input has a 
% duration < Lconv, Matlab performs zero-padding 
% filter_f = fft(h, Lconv);
filter_f = fft([h, zeros(1, Lconv- length(h))]);

for b = 1:num_blocks
    
    block_f = fft([x_pad(1 + (b-1)*(L):(b-1)*L + L), zeros(1, Lconv-L)]);
    
    y_block = ifft(block_f .* filter_f);
    y_oa(1 + (b-1)*(L):(b-1)*L + Lconv) = ...
        y_oa(1 + (b-1)*(L):(b-1)*L + Lconv) + y_block;   
    
end

y_oa = y_oa(1:length(x) + length(h) -1);
hold on; 
stem(n_y, y_oa, '--');

%% overlap and save method.

L = 6;
% number of wrong samples = overlap size. 
error = length(h) - 1;

% pad arrays
h_L = padarray(h, [0, L-length(h)], 'post');
x_L = padarray(x, [0, error], 'pre');

% how many overlapped blocks of length(L) in x_L?
num_blocks = ceil((length(x_L) - L)/(L-error)) + 1;

% in case x_L is not long enough, add zeros
if length(x_L) < (num_blocks - 1)*(L-error) + L
    x_L = [x_L, zeros(1, (num_blocks - 1)*(L-error) + L - length(x_L))];
end

y_os = [];
filter_f = fft(h_L);

for b = 1: num_blocks    
    
    block_f = fft(x_L(1 + (b-1)*(L - error):(b-1)*(L - error) + L));
    y_block = ifft(block_f .* filter_f);
    y_os = [y_os, y_block(length(h):end)];
    
end

y_os = y_os(1:length(x) + length(h) -1 );

hold on, 
stem(n_y, y_os, '-.o')

