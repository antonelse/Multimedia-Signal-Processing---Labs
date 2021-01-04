clear all
close all
clc

%Given the audiofile «flute.wav»
%% 1. Load the signal x
[x, Fs]=audioread('flute.wav');

%% 2. find the pitch frequency of the note ? (in Hertz) by means of the autocorrelation
%       use the function xcorr
[r, lags]=xcorr(x);
r=r(lags>=0);
lags=lags(lags>=0);
first_neg=find(r<0,1,'first');
r(1:first_neg)=-1;
[maxR, maxI]=max(r);
maxL=lags(maxI);
p=Fs/maxL;
fprintf('The pitch frequency in Hertz is %.2f\n',p);

%% 2. Compute the spectrogram of the signal

Nfft=2048;
win=hann(Nfft);
hopsize=Nfft/2;
X_stft=spectrogram(x,win, hopsize);


%% 3. Plot the magnitude of the STFT  for a generic frame ?with the frequency axis in Hertz
generic_frame=42;
f=linspace(0, Fs/2, Nfft/2+1);
figure;
plot(f, abs(X_stft(:,generic_frame)));
xlabel('f [Hz]')
ylabel('20 log |X_{stft}(:, frame)|');
title('Magnitude of the STFT for a generic frame');

%% 4. Find and stem, over the previous figure, the sample corresponding (or closest to) the p
hold on;
%p_i: (Nfft/2+1)= p:Fs/2 
% p_i= p* (Nfft/2+1)/(Fs/2);
p_i = round(p* (Nfft/2+1)/(Fs/2));

stem(f(p_i), abs(X_stft(p_i,generic_frame)));
legend('Magnitude','Pitch sample');