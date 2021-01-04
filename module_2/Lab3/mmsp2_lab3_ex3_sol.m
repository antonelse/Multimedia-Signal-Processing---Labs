clearvars
close all
clc

%% Load the audio file
%!!!!!!! ns deriva da una quantizza del file gb del file precedente

[x, Fs] = audioread('ns.wav');


%% Quantize it using PCM with 1:8 bit, using floor, ceil and round in the quantizer.
%% Compute the SNR.
R = 1:8;

max_x = max(x);
min_x = min(x);

MSE_f = zeros(length(R),1);
MSE_c = zeros(length(R),1);
MSE_r = zeros(length(R),1);

for ii = 1:length(R)

    % Determine delta
    delta = (max_x-min_x)/(2^R(ii)-1);

    % Quantize x, comparando le diverse approssimazioni
    x_f = delta * floor(x/delta) + delta/2; %appross difetto
    x_c = delta * ceil(x/delta) + delta/2; %appross per accesso
    x_r = delta * round(x/delta) + delta/2; %approssimazione

    % Compute MSE_pcm
    MSE_f(ii) = mean((x-x_f).^2);
    MSE_c(ii) = mean((x-x_c).^2);
    MSE_r(ii) = mean((x-x_r).^2);

end
SNR_f = pow2db(var(x)./MSE_f);
SNR_c = pow2db(var(x)./MSE_c);
SNR_r = pow2db(var(x)./MSE_r);

%% Plot the SNR for the different method. What's wrong?
figure
plot(R,SNR_f, R, SNR_c, R, SNR_r);
legend('Floor','Ceil','Round');
grid on;
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');

%%il grafico noto che floor ed il round sono migliori soprattutto per il picco
%dell'SNR a 5 bit che il ceil non rileva. Questo perchè ns era già stato codificato
%con 5bit/sample nel file precedente
