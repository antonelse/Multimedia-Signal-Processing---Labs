clear all
close all
clc
%Given the 16kHz-sampled audiofile Toms_diner_16.wav
% 1. Load the file in the signal x 

[x, Fs]=audioread('Toms_diner_16.wav');
x=x.';

% 2. Design a FIR lowpass filter b with cut-off frequency=6,000 Hz
%       2.1 use the function fir1 and M=32 weights


freq_cutoff=6000;
norm_freq_cutoff=freq_cutoff/(Fs/2);
M=32;
b=fir1(M-1,norm_freq_cutoff);

% 3. Process the signal x with the filter h=b using the Overlap and Add technique
    % 3.1 use block size N=400 choose the kind parameters of the window you prefer
    
h=b;
N=400;
Nfft=2^ceil(log2(N+M-1));
win=bartlett(N).';
hopsize=N/2;
x_=[zeros(1,N/2), x, zeros(1,N/2)];
H=fft(h,Nfft);
k=1;
y=zeros(size(x_));
while k-1+N<length(x_)
   x_f=x_(k:k-1+N).*win;
   X_f=fft(x_f, Nfft);
   y_f=ifft(X_f.*H, Nfft);   
   y(k:k-1+N)=y(k:k-1+N)+y_f(1:N);  
   %% THIS IS WRONG!
   % we are interested in the first
   % N+M-1 samples, we can keep the first
   % Nfft samples, we can NOT only consider
   % the first N sampes
   % please look at the ex8OLALP solution   
   k=k+hopsize-1;
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


