%% Filtering in time and frequency
clear all
clc
close all
% Given 
% the filter h=[-1,0,1] 
% the signal x=[3 1 0 2 1 5] 
% time indices starting both from zero:

h=[-1,0,1];
x=[3 1 0 2 1 5];

%% 1) compute y as the linear convolution of h and x in the time domain
%    * you can use the conv function

y=conv(x,h);

%% 2) compute y2 as the linear convolution of h and x in the frequency domain
%    * you can use the fft function

Nfft=2^ceil(log2(length(x)+length(h)-1));
X=fft(x, Nfft);
H=fft(h, Nfft);
Y=X.*H;
y2=ifft(Y,Nfft);
y2=y2(1:length(x)+length(h)-1);
%% 3) Plot the squared error between y and y2 sample by sample

n=0:length(y)-1;
figure;
plot(n,(y-y2).^2);
xlabel('n');
ylabel('(y-y_2)^2');
title('Squared error sample by sample');
%% 4) Plot the two signals y and y2 together
figure;
stem(n,y);
hold on;
stem(n+0.01,y2,'--');
hold off;
xlabel('n');
ylabel('y, y_2');
title('Output of the linear convolution');

