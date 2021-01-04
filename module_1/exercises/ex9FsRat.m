clear all
close all
clc

% Given the 16kHz-sampled audiofile Toms_diner_16.wav
% 1. Load the file in the signal x

[x, Fs]=audioread('Toms_diner_16.wav');
x=x.';

% 2. Create y from x with frequency rate =12 kHz 
q=12000/Fs;
[L, M]=rat(q);

y_int=zeros(1,length(x)*L);
y_int(1:L:end)=x;
y_int=filter(L*fir1(31,1/L),1,y_int);
y_dec=filter(fir1(31,1/M),1,y_int);
% this is equivalent to
% y_dec=filter(L*fir1(31,min(1/L,1/M)),1,y_int);

y_dec=y_dec(1:M:end);
% 3. Plot two subplots with the original signal x and the modified version 
t_x=[0:length(x)-1]/Fs;
t_y=[0:length(y_dec)-1]/(Fs*q);
figure;
subplot(2,1,1);
plot(t_x,x);
xlabel('time');
ylabel('x');
title('original signal at 16 kHz');
subplot(2,1,2);
plot(t_y,y_dec);
xlabel('time');
ylabel('y');
title('re-sampled signal at 12 kHz');


