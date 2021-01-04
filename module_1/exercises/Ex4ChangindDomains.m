clear all
close all
clc

%Given a filter with H(z)
% z=[0.8, ?0.8, 0.8e^(?j\pi/3),0.8e^(j\pi/3)]
% p=[0+0.5j, 0?0.5j, ?0.6+0.5j, ?0.6?0.5j]

z=[0.8, -0.8, 0.8*exp(-1i*pi/3),0.8*exp(1i*pi/3)].';
p=[0+0.5*1i, 0-0.5*1i, -0.6+0.5*1i, -0.6-0.5*1i].';

%% 1. Plot the z-plane

figure;
zplane(z,p);
title('z-plane');
%% 2. Find the corresponding difference equations B(z) and A(z)

b=poly(z);
a=poly(p);

%% 3. Plot the transfer function and the impulse response of H(z) 
% in 1024 points
N=512;
[H,w]=freqz(b,a,N);
h=filter(b,a,[1,zeros(1,N-1)]);

figure; 
subplot(2,1,1);
plot(w,abs(H));
xlabel('\omega');
ylabel('|H|');
title('Magnitude of the TF');
subplot(2,1,2);
plot(w,angle(H));
xlabel('\omega');
ylabel('\angle H');
title('Phase of the TF');

figure;
plot([0:N-1],h);
xlabel('n');
ylabel('h(n)');
title('impulse response');