clearvars
close all
clc

%% 1) Load the image 'lena512color.tiff'

im = double(imread('lena512color.tiff'));
%faccio double

%% 2) Let x be the red channel and y the green channel of the image.
%% Quantize y with PCM and DPCM (R=1,2,...,8 bits) using:
%%      2.1) y_hat(n) = a*x(n)+b
%%      2.2) y_hat(n) = randn*x(n) + randn*100

x = im(:,:, 1); %red channel (:,:,1)
y = im(:,:, 2); %green channel

x = x(:); % la x è una matrice e facendo così prendo tutti i valori e li "stendo" su un vettore
y = y(:);

% Compute coefficients a and b with LS
coeff = [x ones(size(x,1),1)]\y; %nb "\" risolve l'eq lineare  ax+b=0; gli uni sono riferiti
%al temine noto; nb è una matrice a 2 colonne

%ALTERNATIVA by fotis
%x_ones =[X, ones(length(X),1)];
%coeff = x_ones \ Y;

a = coeff(1); %coeff avrà 2 colonne, termine riferito alla x e al termine noto
b = coeff(2);

% DPCM predictor

% please notice that for the sake of the exercise we are assuming the
% decoder knows x, a and b. In practice this is not a good idea since the
% codec depends on the signal. A more realistic approach would have been
% encoding x_tilde = q(x) with PCM and use x_tilde for the prediction.

y_hat = a*x + b; %y_hat è dato, quindi FISSO e nn si aggiorna ad ogni iterazione !
%Infatti nn ha cicli for e lo mette fuori
%vedi righe dentro il for loop, dove diversamente da quanto facevamo prima, non ho una
%seconda iterazione ad ogni ciclo

% Dummy DPCM predictor
rng(21);
%rng serve per dare la stessa seq randomica ogni volta che invoco randn;
%questo mi serve per eccitare il sist sempre con la stessa seq; comodo per confronti tra metodi
a_d = randn;
b_d = randn*100;
y_hat_d = a_d*x+b_d*100; %il *100 è dato dalla consegna !

R = 1:8;

MSE_dpcm = zeros(length(R), 1);
MSE_dpcm_d = zeros(length(R), 1);
MSE_pcm = zeros(length(R), 1);
for ii = 1:length(R)

    % DPCM
    d = y - y_hat;
    delta_dpcm = (max(d) - min(d))/ (2^R(ii)-1);
    d_tilde = delta_dpcm * floor(d/delta_dpcm) + delta_dpcm/2;
    y_tilde_dpcm = d_tilde + y_hat;

    % PCM
    delta_pcm = (max(y) - min(y))/(2^R(ii)-1);
    y_tilde_pcm = delta_pcm * floor(y/delta_pcm) + delta_pcm/2;

    % Dummy DPCM
    d = y - y_hat_d;
    delta_dpcm = (max(d) - min(d))/ (2^R(ii)-1);
    d_tilde = delta_dpcm * floor(d/delta_dpcm) + delta_dpcm/2;
    y_tilde_dpcm_d = d_tilde + y_hat_d;

    % MSE
    MSE_dpcm(ii) = mean((y - y_tilde_dpcm).^2);
    MSE_pcm(ii) = mean((y - y_tilde_pcm).^2);
    MSE_dpcm_d(ii) = mean((y - y_tilde_dpcm_d).^2);


end

% SNR
SNR_dpcm = pow2db(var(y) ./ MSE_dpcm);
SNR_pcm = pow2db(var(y) ./ MSE_pcm);
SNR_dpcm_d = pow2db(var(y) ./ MSE_dpcm_d);

%% Compare the R-D curves
figure
plot(R, SNR_dpcm, R, SNR_pcm, R, SNR_dpcm_d);
legend('DPCM', 'PCM', 'Dummy DPCM')
grid on;
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');
