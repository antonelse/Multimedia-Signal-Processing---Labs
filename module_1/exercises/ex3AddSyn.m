%% ex3 Addythive synthesis
% Open the wav file flute.wav as ? and try to synthetically (and empirically) replicate it as a sum of K=3 sinusoids

clear all
close all
clc

[x, Fs]=audioread('flute.wav');
%% 1. Compute and plot the spectrogram of x

Nfft=2048;
R=Nfft/2;
f=linspace(0,Fs/2,Nfft/2+1);
X=spectrogram(x,bartlett(Nfft),R,Nfft);
t=linspace(0,length(x)/Fs,size(X,2));

figure; imagesc(t,f,20*log10(abs(X)));
xlabel('Time [s]');
ylabel('Frequency [Hz]')
%% 2. From a random frame of the spectrogram, find the K main sinusoids 
% using the provided function findPeaks

frame=32;
K=3;
peaks=findPeaks(abs(X(:,frame)),K);
figure; plot(f,abs(X(:,frame)));
hold on; stem(f(peaks), abs(X(peaks,frame))); hold off

%% 3. Find the frequencies f_k and, from the STFT, 
% the (time-variant) magnitudes A_k(t)
ampl_STFT=abs(X(peaks,:));
ampl_STFT=ampl_STFT/max(ampl_STFT(:));
figure; plot(ampl_STFT');



%% 4. Compute y=sum_k A_k*cos(2*pi*f_k*t)

y=zeros(size(x));
freqs=f(peaks)*2*pi;
n=0:length(x)-1;
n=(n')/Fs;

for k=1:length(peaks)
   A_k=interp([0,ampl_STFT(k,:),0], R);
   % Interpolate A_k using the function interp 
   % and L=hopsize of the spectrogram window
   
   A_k=A_k(1:length(n));
   A_k=A_k.';   
   y_k=A_k.*cos(freqs(k)*n);
   y=y+y_k;
   
end

%% 5. Normalize it and listen to the result

y=y/max(abs(y));
sound(x,Fs);
pause(length(x)/Fs+1);
sound(y,Fs);

