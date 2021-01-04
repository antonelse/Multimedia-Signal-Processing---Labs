%% MMSP2 - Lab 3
%  Exercise 2 - Predictive coding

clear
close all
clc

%% Load the stereo file ‘mso.wav’ and define xl and xr as
%% the left and right channels, respectively.

[x, Fs] = audioread('mso.wav');

N = size(x, 1);
xl = x(:, 1);
xr = x(:, 2);
% x sono 2 colonne poichè segnale stereo

%% 1) Build a DPCM codec. Use the left channel xl as a signal and
%%      1.1) xl(n-1)
%%      1.2) xr(n)
%%      1.3) dummy 5*xl(n)
%% as predictor.
%%    Use PCM to initialize the codec.
R = 1:8;

MSE_l = zeros(length(R),1); %caso 1.1
MSE_r = zeros(length(R),1); %caso 1.2
MSE_d = zeros(length(R),1); %caso 1.3
MSE_pcm = zeros(length(R),1); %che ci serve per il primo sample

for ii = 1:length(R)
    xl_l_tilde = zeros(N,1);
    xl_r_tilde = zeros(N,1);
    xl_d_tilde = zeros(N,1);

    % first sample: PCM
    max_xl = max(xl);
    min_xl = min(xl);
    delta_pcm = (max_xl-min_xl)/(2^R(ii)-1);

    xl_l_tilde(1) = delta_pcm * floor(xl(1)/delta_pcm) + delta_pcm/2;

    % next samples: DPCM
    d_l_aux = xl(2:end) - xl(1:end-1);
    % without quantizer, this would be the residual,
    %poichè il predittore in questo caso è se stesso ma n-1 quindi parte da 1 fino a N-1 (cioè end-1)
    %PRIMA (con AR process) lo assumevo come rumore. Ora nn posso.
    %NOTA BENE. La definizione di d dipende dal predictor !!! SOLO QUESTA PARTE CAMBIA (aka la parte
    %di definizione di d tilde di ogni casistica è in base al relativo x_hat)
    delta_l = (max(d_l_aux) - min(d_l_aux))/(2^R(ii)-1);

    d_r_aux = xl - xr;
    delta_r = (max(d_r_aux) - min(d_r_aux))/(2^R(ii)-1);

    d_d_aux = xl - 5*xl;
    delta_d = (max(d_d_aux) - min(d_d_aux))/(2^R(ii)-1);

    for nn = 2:N
        % predict xl(n) from xl(n-1)
        x_l_hat = xl_l_tilde(nn - 1);
        d_l = xl(nn) - x_l_hat;
        d_l_tilde = delta_l * floor(d_l/delta_l) + delta_l/2;
        xl_l_tilde(nn) = d_l_tilde + x_l_hat;

        % predict xl(n) from xr(n)
        x_r_hat = xr(nn);
        d_r = xl(nn) - x_r_hat;
        d_r_tilde = delta_r * floor(d_r/delta_r) + delta_r/2;
        xl_r_tilde(nn) = d_r_tilde + x_r_hat;

        % predict xl(n) from dummy linear combination
        x_d_hat = 5*xl_l_tilde(nn);
        d_d = xl(nn) - x_d_hat;
        d_d_tilde = delta_d * floor(d_d/delta_d) + delta_d/2;
        xl_d_tilde(nn) = d_d_tilde + x_d_hat;
    end

    % Let's build a PCM codec, just for comparison, CIOE' NO LOOP, e nemmeno l'open loop
    xl_pcm_tilde = delta_pcm * floor(xl/delta_pcm) + delta_pcm/2;

    % MSE
    MSE_l(ii) = mean((xl-xl_l_tilde).^2);
    MSE_r(ii) = mean((xl-xl_r_tilde).^2);
    MSE_d(ii) = mean((xl-xl_d_tilde).^2);
    MSE_pcm(ii) = mean((xl - xl_pcm_tilde).^2);


end

SNR_l = pow2db(var(xl)./MSE_l);
SNR_r = pow2db(var(xl)./MSE_r);
SNR_d = pow2db(var(xl)./MSE_d);
SNR_pcm = pow2db(var(xl)./MSE_pcm);

%% 2) Compare R-D curves
figure
plot(R,SNR_l,R,SNR_r,R,SNR_d,R, SNR_pcm);
legend('$\hat{x} = xl(n-1)$','$\hat{x} = xr(n)$','$\hat{x} = 5*xl(n)$', 'PCM', 'Interpreter', 'latex');
%$\ è per fare il carattere speciale hat
grid on;
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');
