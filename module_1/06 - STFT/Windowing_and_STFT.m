clear all
close all
clc

%% Windows: first try

f=510; N=512; Fs=10000;
n=(0:N-1)/Fs;
x=exp(1i*2*pi*f*n);
x=x.';
w=window(@bartlett,50);
w(N)=0;
xw=x.*w;

X=fft(x,N);
W=fft(w,N);
XW=fft(xw,N);
F=linspace(0,Fs,N);

figure;
subplot(3,1,1)
plot(n,real(x));
xlabel('n'); ylabel('x(n)'); title('Signal');
subplot(3,1,2)
plot(n,w);
xlabel('n'); ylabel('w(n)'); title('Window');
subplot(3,1,3)
plot(n,real(xw));
xlabel('n'); ylabel('x_w(n)'); title('Windowed signal');


figure;
subplot(3,1,1)
plot(F,abs(X));
xlabel('frequency [Hz]'); ylabel('|X(k)|'); title('Signal');
subplot(3,1,2)
plot(F,abs(W));
xlabel('frequency [Hz]'); ylabel('|W(k)|'); title('Window');
subplot(3,1,3)
plot(F,abs(XW));
xlabel('frequency [Hz]'); ylabel('|X_w(k)|'); title('Windowed signal');

%% Causal windows
N=101;
win=window(@bartlett,N);
win=[zeros(N,1); win; zeros(N,1)];
n=linspace(-N/2,N/2,N*3);
win_shifted=win(N+1:end);
win_shifted(3*N)=0;


[WIN, w]=dft(win,n); 
% need to use a modified DFT to see with negative n
WIN_shifted=dft(win_shifted, n-min(n));


figure; 
subplot(2,2,1);
plot(n, win);
ylabel('w(n)'); xlabel('n');
title('non-causal window');

subplot(2,2,2);
plot(n-min(n), win_shifted);
ylabel('w(n)'); xlabel('n');
title('causal window');

subplot(2,2,3);
plot(w/pi,angle(WIN));
ylim([-2*pi,2*pi]);
ylabel('\angle W(z)'); xlabel('normalized \omega');
title('zero-phase');

subplot(2,2,4);
plot(w/pi,unwrap(angle(WIN_shifted)));
%ylim([-2*pi,2*pi]);
ylabel('\angle W(z)'); xlabel('normalized \omega');
title('linear-phase');


%% Types of windows
N_win=64;
N=1024;

rect_win=window(@rectwin,N_win);
show_win(rect_win, N, 'rect');

triang_win=window(@bartlett,N_win);
show_win(triang_win, N, 'triang');

hann_win=window(@hann,N_win);
show_win(hann_win, N, 'Hann');

hamming_win=window(@hamming, N_win);
show_win(hamming_win, N, 'Hamming');

black_win=window(@blackman, N_win);
show_win(black_win, N, 'Blackman');

kaiser_win=window(@kaiser, N_win);
show_win(kaiser_win, N, 'Kaiser');


gauss_win=window(@gausswin, N_win);
show_win(gauss_win, N, 'Gaussian');


%%
clear all

N=2^14;
n=(0:N-1).';
omega1=2*pi*1000/N;
omega2=2*pi*1000/N + 2*pi/40;
DeltaOmega=abs(omega1-omega2);
x=cos(omega1*n)+cos(omega2*n);

omega_norm=linspace(-1,1,N);    
M=[16, 32, 64, 128];
figure;
for m =1:length(M)
    subplot(2,2,m);
    M_i=M(m);
    w_R=window(@rectwin,M_i)/M_i;
    w_R(N)=0;
    y=x.*w_R;
    Y=fft(y,N);
    Y=fftshift(Y);
    plot(omega_norm, 20*log10(abs(Y)));
   
   % plot(omega_norm, 20*log10(abs(Y)));
    xlabel('\omega [\pi units]');
    ylabel('|Y(k)|');
    title(['M= ' num2str(M_i)]);
end

%% Minimum M
clear all
N=2^14; n=(0:N-1).';
w_norm=linspace(-1,1,N);    
omega1=0.2*pi; omega2=0.24*pi;
x=exp(1i*omega1*n)+exp(1i*omega2*n);
Delta=abs(omega1-omega2);
L=2; M=2*L*2*pi/Delta; % rectWin
M=ceil(1.1*M); % adding a 10% of the value
w_R=window(@rectwin,M)/M; w_R(N)=0; 
y=x.*w_R; Y=fftshift(fft(y,N));
figure;
plot(w_norm, 20*log10(abs(Y)));
xlabel('\omega [\pi units]');
ylabel('|Y(k)|');
title('Equal amplitude signals');

%% Different amplitude signals
x2=exp(1i*omega1*n)+0.1*exp(1i*omega2*n);
%w_R=window(@hann,2*M); w_R(N)=0; 

y2=x2.*w_R; Y2=fftshift(fft(y2,N));
figure;
plot(w_norm, 20*log10(abs(Y2)));
xlabel('\omega [\pi units]');
ylabel('|Y(k)|');
title('Different amplitude signals');


%% Different amplitude signals with Hamming window
L=2; M=2*L*2*pi/Delta; % rectWin
M=ceil(1.1*M); % adding a 10% of the value
w_H=window(@hamming,M)/M; w_H(N)=0; 
y3=x2.*w_H;

Y3=fftshift(fft(y3,N));
figure;
plot(w_norm, 20*log10(abs(Y3)));
xlabel('\omega [\pi units]');
ylabel('|Y(k)|');
title('Different amplitude signals with Hamming window');
 
%% overlap and add
clear all
[x, Fs]=audioread('gb.wav');
N=length(x);
K=1000;
delta=[1;0];
delta(K)=0;
h=filter(1,[1 -0.99],delta);

M=floor(0.050*Fs); % frame of 50ms
w=window(@bartlett,M);
R=floor(M*0.5); % 50% overlap

Nfft=2^ceil(log2(K+M-1));

y_=conv(x,h);

h(Nfft)=0;
H=fft(h, Nfft);

%% is COLA respected?
times=10;
COLA=zeros((M-R)*times+R,1);
w=window(@rectwin,M);
figure;
hold on
plot(COLA);
pause(0.5);
for i=0:times-1
   COLA_win=zeros(size(COLA));
   COLA((i*R)+1:(i*R)+M)=COLA((i*R)+1:(i*R)+M)+w; 
   COLA_win((i*R)+1:(i*R)+M)=w;   
   plot(COLA_win);  
   plot(COLA,'k');
   ylim([0,5]);
   pause(1.5);
end

clear COLA COLA_win
%%

x(N+M+R+Nfft)=0;
y=zeros(size(x));
%x(N+M+R)=0;
for m=1:floor(N/R)+1
    mR=(m-1)*R+1;
    x_m=x(mR:mR+M-1);
    x_wm=x_m.*w;
    x_wm(Nfft)=0;
    X_wm=fft(x_wm);
    Y_m=X_wm.*H;
    y_m=ifft(Y_m);
    y(mR:mR+Nfft-1)=y(mR:mR+Nfft-1)+y_m;   
end

y=y(1:N+K-1);
figure; plot(abs(y-y_));


%% Showing OLA
clear all
[x, Fs]=audioread('gb.wav');
N=length(x);
K=1000;
delta=[1;0];
delta(K)=0;
h=filter(1,[1 -0.99],delta);

M=floor(0.050*Fs); % frame of 50ms
w=window(@bartlett,M);
R=floor(M*0.5); % 50% overlap

Nfft=K+M-1;

y_=conv(x,h);

h(Nfft)=0;
H=fft(h, Nfft);
x=x(R*15:R*25);
N=length(x);
x(N+M+R)=0;
y=zeros(size(x));

figure;
subplot(4,2,2);
plot(0:K-1,h(1:K));
xlabel('n');
ylabel('h(n)');
title('h(n)');
subplot(4,2,4);
plot(abs(H));
xlabel('k');
ylabel('|H(k)|');
title('|H|');

for m=1:floor(N/R)+1
    subplot(4,2,1);
    plot(x);
    title('x')
    mR=(m-1)*R+1;
    x_m=x(mR:mR+M-1);
    hold on
    x_=zeros(size(x));
    x_(mR:mR+M-1)=w;
    plot(x_,'r');
    hold off
    
    %subplot(5,2,3);
    %x_=zeros(size(x));
    %x_(mR:mR+M-1)=x_m;
    %plot(x_);
    %title('x_m')
    x_wm=x_m.*w;
    
    subplot(4,2,3);    
    x_(mR:mR+M-1)=x_wm;
    plot(x_);
    title('x_{w,m}');
    
    
    x_wm(Nfft)=0;
    X_wm=fft(x_wm);
    
    subplot(4,2,6);    
    plot(abs(X_wm));
    title('|X_{w,m}|');
    
    Y_m=X_wm.*H;
    
    subplot(4,2,8);    
    plot(abs(Y_m));
    title('|Y_m|=|X_{w,m}| |H|');
    
    y_m=ifft(Y_m);
    y(mR:mR+Nfft-1)=y(mR:mR+Nfft-1)+y_m;   
    
    subplot(4,2,5);
    y_=zeros(size(y));
    y_(mR:mR+Nfft-1)=y_m;
    plot(y_);
    title('y_m');
    subplot(4,2,7);    
    plot(y);
    title('y');
    if m==1
        pause(); 
    else
        pause(3);
    end
end

%% STFT with audio
clear all
[x, Fs]=audioread('gb.wav'); 
N=length(x);
M=floor(0.050*Fs); 
R=floor(M*0.5);
w=window(@bartlett,M);
Nfft=M; 
x(N+M)=0;
X=zeros(floor(Nfft/2+1),floor(N/R)+1);
for m=1:floor(N/R)+1
  mR=(m-1)*R+1; x_m=x(mR:mR+M-1);    
  x_wm=x_m.*w; X_wm=fft(x_wm, Nfft);
  X(:,m)=X_wm(1:floor(Nfft/2+1));
end
X=flipud(X);
time=(1:size(X,2))*0.025;
freq=linspace(Fs/2,0,size(X,1));

figure; imagesc(time, freq,20*log10(abs(X)));
title('STFT X [dB]')
xlabel('Time');
ylabel('Freq');
axis xy

%% A changing frequency signal
close all
n=[0:2^14];
x=cos(2*pi*n/100);
x(n>2^13)=cos(2*pi*n(n>2^13)/200);

figure;
plot(n,x);
figure;
plot(abs(fft(x,2^14)));

figure;
Nspec=1024;
imagesc(abs(spectrogram(x,window(@bartlett, Nspec),Nspec/2,Nspec)));


%% invert magnitude and phase in spectrograms
[y_fgi,Fs]=audioread('fgi.wav');
[y_wttj,Fs]=audioread('wttj.wav');
y_fgi=mean(y_fgi.');
y_wttj=mean(y_wttj.');
%mean converts mono to stereo
logN=floor(log2(0.05*Fs));
N=2^logN;
Noverlap=2^(logN-1);

Nmin=min(length(y_fgi),length(y_wttj));
y_fgi=y_fgi(1:Nmin);
y_wttj=y_wttj(1:Nmin);

win=window(@hamming,N);
Y_fgi=spectrogram(y_fgi, win,Noverlap);
Y_wttj=spectrogram(y_wttj, win,Noverlap);
Y_fm_wa=abs(Y_fgi).*exp(1i*angle(Y_wttj));
Y_wm_fa=abs(Y_wttj).*exp(1i*angle(Y_fgi));
y_fm_wa=istft(Y_fm_wa,Noverlap);
y_wm_fa=istft(Y_wm_fa,Noverlap);

sound(y_fm_wa,Fs)
sound(y_wm_fa,Fs)
