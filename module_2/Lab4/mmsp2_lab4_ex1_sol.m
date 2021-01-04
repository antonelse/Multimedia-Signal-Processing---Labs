%% MMSP2 - Lab 4
%  Exercise 1 - Transform coding

clear
close all
clc

%% 1) Load the first 4s of the file 'gb.wav' and quantize it with PCM and R=8 bit.
%%    Compute the MSE and perceptually evaluate the result.
[x, Fs] = audioread('gb.wav');

R = 8;
len = 4;
x = x(1:len*Fs);

delta_pcm = (max(x)-min(x))/(2^R-1);
x_pcm = delta_pcm*floor(x/delta_pcm) + delta_pcm/2;

mse_pcm = var(x-x_pcm);
%analogo di mse_pcm = mean((x - x_pcm)^2);
snr_pcm = pow2db(var(x)/mse_pcm);

disp(['SNR PCM: ' num2str(snr_pcm) 'dB' ]); %disp lo fa spuntare nella command wind
% Listen the audio track
sound(x, Fs)
sound(x_pcm, Fs)


%% 2) Consider groups of 8 symbols and quantize them using an optimal allocation of the 8 bits
% complete the function transform_coding
% hint: consider each block of 8 symbols as if composed by 8 transform
% coefficients

N = 8; %group of symbols length

% Define the transformation matrix
T_eye = eye(N); %eye: matrice identità, cioè uno sulla diag
[R_eye,snr_eye] = transform_coding(x,T_eye,R); %questa funz riceve in ingresso il segnale, la matr di trasf e il num dei bit.
%Restituisce l'allocazione perfetta dei bit e il snr associato

disp(['SNR Eye: ' num2str(snr_eye) 'dB' ]);

%% 3) Consider DCT transformation and repeat step 2 over transformed
%%    coefficients. Find the distortion and evaluate the perceived quality.

N = 8; %group of symbols length

% Define the transformation matrix, vedi l'hint !!! è definita così
T_dct = zeros(N,N);
T_dct(1,:) = sqrt(1/N);

l = 1:N;
for k = 2:N
    T_dct(k,:) = sqrt(2/N)*cos(pi/(2*N)*(k-1)*(2.*l-1));
end

%see also: T = dctmtx(N);

[R_dct,snr_dct] = transform_coding(x,T_dct,R);
%quindi invece di applicare la matr di trasf unitaria (T), applico la matrice di trasf della DCT

disp(['SNR DCT: ' num2str(snr_dct) 'dB' ]);

%% 3) Consider a Karhunen-Loeve transformation and repeat step 2 over transformed
%%    coefficients. Find the distortion and evaluate the perceived quality.

% compute correlation !!!!!!!
% hint: remove the average from the signal
% hint: compute many 8x8 correlation matrix, then average them

% remove signal mean
x_zm = x - mean(x);

% Estimate the autocorrelation function for groups of 8 symbols
X_zm = reshape(x_zm,N,length(x_zm)/N);
%reshape(X,M,N) or reshape(X,[M,N]) returns the M-by-N matrix whose elements are taken columnwise from X.
%quindi credo matrice dal segnale vett colonna. Ogni colonna ha N elementi e rappresenta un chunks di segnale !

RR = zeros(N,N,length(x_zm)/N);
%creo un array in 3D, così posso mediare le varie matrici di correlazione 2D che si creeranno autocorrelando i chunks
for ii = 1:length(x_zm)/N %aka 3rd dimension, cioè io numero di matrici di corr che medierò insieme; aka layer
    RR(:,:,ii) = X_zm(:,ii)*X_zm(:,ii)';
    %non faccio il mean (aka val atteso) perchè lo faccio direttamente alla fine tra tutte le matr di corr
end
RR_mean = mean(RR,3);
%medio le matrici 2D, quindi le prendo tutte (che sono lungo la 3rd dimension). R_mean è una 2D

% compute eigenvectors and eigenvalues of the autocorrelation matrix
[V,~] = eig(RR_mean);
%returns diagonal matrix D (aka ~ in questo caso, poichè non è di nostro interesse) of eigenvalues and matrix V
%whose columns are the corresponding right eigenvectors

% define the transformation matrix (BE AWARE THAT WE NEED A TRANSFORMATION THAT THE ROWS
%(e non colonne, ovvero come ce li dà il comando eig() !! ) ARE OUR PROJECTION BASIS !!!), quindi faccio la trasposizione
T_klt = V';

[R_klt,snr_klt] = transform_coding(x_zm,T_klt,R); % of course if we want to go back to x we need to add the mean back

disp(['SNR KLT: ' num2str(snr_klt) 'dB' ]);


%% 4) Plot the amount of bits allocated for each of the proposed solutions

figure();
subplot(3,1,1);
bar(R_eye); %grafico a barre
xlabel('Coefficient');
ylabel('bit');
title('Eye');
grid on;
%la trasf unitaria infatti alloca gli stessi bit per ogni coeff

subplot(3,1,2);
bar(R_dct);
xlabel('Coefficient');
ylabel('bit');
title('DCT');
grid on;

subplot(3,1,3);
bar(R_klt);
xlabel('Coefficient');
ylabel('bit');
title('KLT');
grid on;


IMPORTANTE

%la basis function (righe di T) moltiplicate per le colonne di x (che rappresentano ognuna
%un gruppo di simboli lungo N), vengono rappresentate da un singolo coefficiente della matrice Y.
%La prima colonna della matrice Y mi indica quanto il primo gruppo sia CL di tutte le
%basis function (righe) della mia matrice di trasf T
oppure
%la prima riga di Y i indica quanta influenza della prima base function c'è stata nella composizione
%di ogni gruppo di simboli
