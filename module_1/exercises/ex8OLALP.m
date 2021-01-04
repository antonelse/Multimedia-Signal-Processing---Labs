clear all
close all
clc

%Given the 16kHz-sampled audiofile Toms_diner_16.wav
%% 1. Load the file in the signal x
[x, Fs]=audioread('Toms_diner_16.wav');
%% 2. Design a FIR lowpass filter b 
%           with cut-off frequency=6,000 Hz
% 2.1. use the function fir1 and M=32 weights
M=32;
b=fir1(M-1,6000/(Fs/2));

%% 3. Process the signal x with the filter h=b %
% using the Overlap and Add technique
% 3. 1. use block size N=400 choose the parameters and 
% kind of the window you prefer
N=400; h=b.';
Nfft=2^ceil(log2(N+M-1));
win=window(@bartlett, N); hopsize=N/2;
x_=[zeros(hopsize,1); x; zeros(hopsize,1)];

H=fft(h,Nfft); y=zeros(size(x_));
y=[y; zeros(Nfft,1)];
for mR=1:hopsize:length(x_)-N  
   x_block=x_(mR:mR+N-1).*win;
   X_block=fft(x_block,Nfft);
   Y_block=X_block.*H;
   y_block=ifft(Y_block);
   y(mR:mR+Nfft-1)=y(mR:mR+Nfft-1)+y_block;
end

%% 4. Plot two subplots with the original signal x 
% and the filtered version

figure;
subplot(2,1,1); plot(x);
subplot(2,1,2); plot(y);

