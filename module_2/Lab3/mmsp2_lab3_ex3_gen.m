clearvars
close all
clc

%% Load the file

[x, Fs] = audioread('gb.wav');

%% Quantize x with 5 bit
R = 5;

max_x = max(x);
min_x = min(x);

%prova a fare i casi con diversi ma è solo 1, cioè length(R) è 1
for ii = 1:length(R)

    % Determine delta
    delta = (max_x-min_x)/(2^R(ii)-1);

    % Quantize x
    x_q = delta * floor(x/delta) + delta/2;

end

%% Compare x and x_q
figure
subplot(211);
stft(x);
xlabel('x');
subplot(212);
stft(x_q);
xlabel('x_q');

figure
subplot(211);
histogram(x);
xlabel('x');
subplot(212);
histogram(x_q);
xlabel('x_q');

%% Save the quantized version

audiowrite('ns.wav', x_q, Fs);
