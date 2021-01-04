clear all
close all
clc


%% Polyphase decimation
clear all
[x,Fs]=audioread('Toms_diner_16.wav');
n_x=[0:length(x)-1]/Fs;
M=4;
h=fir1(32,1/M);

e=decompose_filter(h,M);
x_m=decompose_signal(x,M);

y=zeros(size(x_m(:,1)));
for m=1:M
   y_m=filter(e(:,m),1,x_m(:,m));
   y=y+y_m;
end
n_y=[0:length(y)-1]/(Fs/M);

y_dec=filter(h,1,x);

y_dec=y_dec(1:M:end);
n_y_dec=[0:length(y_dec)-1]/(Fs/M);
show_dec=1;
figure;
plot(n_x,x);
hold on;
plot(n_y,abs(y),'r--');

xlabel('t'); ylabel('x(n), y(n)');
if show_dec
    plot(n_y_dec,-abs(y_dec),'g--');
    legend('x(n)','y(n) polyphase', 'y(n) decimated')
else
    legend('x(n)','y(n) polyphase')
end

title('Polyphase decimation');
hold off;
n1=min(length(y),length(y_dec));
figure; 
plot(n_y(1:n1),abs(y(1:n1)-y_dec(1:n1)));
xlabel('time');
ylabel('absolute difference');
title('Difference between normal and polyphase decimation');
%% Polyphase decimation: comparing spectra

% let's choose one frame
Nfft=2048;
%i=ceil(rand(1)*(length(y)-2048)); %random pickup
i=find(n_y>0.58,1,'first'); % look at the bad represantation
%[X,w]=dft_mat(x(floor(i/2):floor(i/2)-1+Nfft*2).', Nfft*2);
[Y,w]=dft_mat(y(i:i-1+Nfft)',Nfft);
[Y_dec,w]=dft_mat(y_dec(i:i-1+Nfft)',Nfft);

figure;
%plot(w/pi,abs(X(1:length(X)/2))/2); hold on;
plot(w/pi,abs(Y)); hold on;
plot(w/pi, abs(Y_dec),'--'); hold off;
xlabel('\omega/\pi'); ylabel('|Y(k)|');
title('Spectrum of polyphase decimation');
legend('polyphase decimation','decimation');

%% Polyphase interpolation
y_int=zeros(length(y)*M+M,1);
for m=1:M
  y_m_f=filter(M*e(:,m),1,y); % scaling
  y_int(m:M:end-M)=y_int(m:M:end-M)+y_m_f;
end

y_int_classic=zeros(length(y)*M,1);
y_int_classic(1:M:end)=y;
y_int_classic=filter(M*h,1,y_int_classic);
fprintf('Difference between interpolation and polyphase interpolation is: %f\n',...
        sum(abs(y_int(1:end-2)-y_int_classic)));
%% Polyphase interpolation: comparing spectra

% let's choose one frame
Nfft=512;
%i=ceil(rand(1)*(length(y)-2048)); %random pickup
i=find(n_x>0.58,1,'first'); % look at the bad represantation
[Y_int,w]=dft_mat(y_int(i:i-1+Nfft)',Nfft);
[Y_int_classic,w]=dft_mat(y_int_classic(i:i-1+Nfft)',Nfft);

figure;
plot(w/pi,20*log10(abs(Y_int))); hold on;
plot(w/pi, 20*log10(abs(Y_int_classic)),'r--'); hold off;
xlabel('\omega/\pi'); ylabel('|Y(k)|');
title('Spectrum of polyphase interpolation');
legend('polyphase interpolation','simple interpolation');

