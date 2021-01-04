clear all
close all
clc

% Given the 16kHz-sampled audiofile Toms_diner_16.wav
% 1. Load the file in the signal x
[x, Fs]=audioread('Toms_diner_16.wav');
x=x.';
% 2. Design a FIR lowpass filter b with cut-off 
%       frequency=6,000 Hz
% 2.1 use the function fir1 and M=32 weights
freq_cutoff=6000;
norm_freq_cutoff=freq_cutoff/(Fs/2);
M=32;
b=fir1(M-1,norm_freq_cutoff);

% 3. Process the signal x with the filter h=b using the Overlap and Save technique
% 3.1 Use block size N=1024

h=b;
N=1024;
x_=[zeros(1,M-1), x, zeros(1,N-1)];
h_=h;
h_(N)=0;
H_=fft(h_,N);
k=1;
y=[];%zeros([]);
while k-1+N<length(x_)
   x_block=x_(k:k-1+N);
   X_block=fft(x_block, N);
   y_block=ifft(X_block.*H_,N);
   y=[y, y_block(M:end)];  
   k=k+N-(M-1);
end
y=y(1:length(x)+length(h)-1);
% 4. Plot two subplots with the original signal x and the filtered version 
t=[0:length(x)-1]/Fs;
figure;
subplot(2,1,1);
plot(t,x);
xlabel('time [s]');
ylabel('x(t)');
title('original signal');
t=[0:length(y)-1]/Fs;
subplot(2,1,2);
plot(t,y);
xlabel('time [s]');
ylabel('y(t)');
title('filtered signal');


