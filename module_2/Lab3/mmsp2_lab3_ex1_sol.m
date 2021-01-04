%% MMSP2 - Lab 3
%  Exercise 1 - Predictive coding

clear
close all
clc


%% 1) Generate 10000 samples of the random process
%%    x(n) = rho*x(n-1) + z(n), where rho=0.95 and z(n)~N(0,0.1)

N = 10000;
rho = 0.95;
z_var = 0.1;
z = randn(N,1) * sqrt(z_var);
x = filter(1,[1 -rho], z); %% la x è l uscita e z è l ingr


%% 2) Build a PCM codec. Quantize the signal with a uniform quantizer and
%%    R bits. Compute the R-D curve for R=1,2,...,8 bits (in terms of SNR)
R = 1:8;

max_x = max(x);
min_x = min(x);

MSE_pcm = zeros(length(R),1);
for ii = 1:length(R)

    % Determine delta
    delta = (max_x-min_x)/(2^R(ii)-1);
    % il -1 nn è a exp, è la stessa cose di dividere per M-1;
    %M era 2^R; è l'alternativa a unique

    % Quantize x
    x_q = delta * floor(x/delta) + delta/2; %MID RISE (repro contiene lo zero, aka ho il gradino in zero)
    %la x sarebbe la solita y, infame

    % Compute MSE_pcm
    MSE_pcm(ii) = mean((x-x_q).^2);

end
SNR_pcm = pow2db(var(x)./MSE_pcm); %MSE_pcm è un vect (quindi anche SNR_pcm), var(x) no


%% 3) Build a predictive codec in open loop. Use the optimal MMSE (Minimum Mean Squared Error) predictor.
%%    Use PCM to initialize the codec
% Remark:
% The optimal MMSE predictor is
% x_hat(n) = rho*x(n-1), hence the prediction residual is
% d(n) = x(n) - x_hat(n) = x(n) - rho*x(n-1) = z(n), aka il rumore ! Ricordati MIDA1
%
% For the first sample of the process (i.e. n = 1) we would have
% d(1) = x(1) - rho*x(0), which is impossible to compute. So, in order to
% highlight this fact, we set d(1) = NaN; please notice that the value of d(1)
% will remain unused.
%
% For the second sample of the process (i.e. n = 2) we have
% d(2) = x(2) - rho*x(1) = z(2).
%
% From that remark follows the definition of the vector d as
d = [NaN; z(2:N)];
%z è stata definita sopra e d è zeta nel caso dell'optimal predictor

max_d = max(d);
min_d = min(d);

MSE_olpc = zeros(length(R),1); %aka dell'OpenLoop
for ii = 1:length(R)
    x_tilde = zeros(N,1);

    % first sample: PCM
    delta_pcm = (max_x-min_x)/(2^R(ii)-1);
    x_tilde(1) = delta_pcm * floor(x(1)/delta_pcm) + delta_pcm/2;

    % next samples: OLPC
    delta_olpc = (max_d-min_d)/(2^R(ii)-1);
    d_tilde = delta_olpc * floor(d/delta_olpc) + delta_olpc/2;

    for nn = 2:N
        x_hat = rho*x_tilde(nn-1)
        x_tilde(nn) = d_tilde(nn) + x_hat; %rho is known also at the receiver!
        % è il predittore, nb rivedi schema sul quad !! e' l'OPEN LOOP !! ha il pred all'inizio
        %(solo primo) e alla parte del ricevitore
    end

    % MSE
    MSE_olpc(ii) = mean((x-x_tilde).^2);
    %l'errore è associato ad un segnale meno il suo quantizzato; ma in questo caso
    %quantizzo la diff tra x e x_hat cioè d, quindi il segnale ricostruito grazie alla
    %quantizz di d e x_tilde quindi uso questo. E' la stesso caso di valutare d- d_tilde (cioè d quantizzato)


end
SNR_olpc = pow2db(var(x)./MSE_olpc);


%% 4) Build a DPCM codec (aka closed loop). Use the optimal MMSE predictor.
%%    Use PCM to initialize the codec.

MSE_dpcm = zeros(length(R),1);
for ii = 1:length(R)
    x_tilde = zeros(N,1);

    % first sample: PCM
    delta_pcm = (max_x-min_x)/(2^R(ii)-1);
    x_tilde(1) = delta_pcm * floor(x(1)/delta_pcm) + delta_pcm/2;
    %nb, fin qui uguale a sopra

    % next samples: DPCM
    delta_dpcm = (max_d-min_d)/(2^R(ii)-1);
    % here we are assuming we know the max and min of d,
    % but in practice we should guess delta_dpcm

    for nn = 2:N %ricorda che questa parte è al decoder !!! che è uguale a quello del PCM
        x_hat = rho*x_tilde(nn-1); %quindi necessito di x_tilde(1) per iniziallizzarlo; questo processo lo faccio fuoti
        d = x(nn) - x_hat; %qui d si aggiorna ogni volta
        d_tilde = delta_dpcm * floor(d/delta_dpcm) + delta_dpcm/2;
        x_tilde(nn) = d_tilde + x_hat;
    end
    % qui uso Closed Loop !! pag 44 tagliasacchi
    % MSE
    MSE_dpcm(ii) = mean((x-x_tilde).^2);
end

SNR_dpcm = pow2db(var(x)./MSE_dpcm);

%% 5) Compare R-D curves for PCM, open-loop DPCM and closed-loop DPCM

plot(R,SNR_dpcm,R,SNR_olpc,R,SNR_pcm); %solo 3 curve, con l asse x sempre R dei bit
legend('DPCM Closed Loop','DPCM Open Loop','PCM');
grid on;
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');
