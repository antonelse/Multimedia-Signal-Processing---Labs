%% MMSP2 - Lab 2
%  Exercise 2 - Uniform and optimal quantizer

clearvars
close all
clc

%% 1) Generate a 10000-sample realization of s_g(n)~N(0,2) and s_u(n)~U with variance 2 and mean 0
%%    hint: use the functions randn() and rand()
var_s_g = 2;
var_s_u = 2;

% Normal distribution
s_g = randn(10000,1) * sqrt(var_s_g);

% Uniform distribution
%%formula nota poiche varianza_unif=(b-a)^2/12, b-a = delta_u
delta_u = sqrt(var_s_u * 12);
%rescaling, aka centratura della distr unif in zero
s_u = rand(10000,1) * delta_u - delta_u/2;

%% 2) Quantize s_g(n) and s_u(n) with M=[4,8,16,32,64,128] levels (cioè 2, 3, 4, 5, 6, 7 bit) and uniform
%%    quantizer. Plot R-D curve for each number of levels. Compare with the
%%    theoretical distortion for a uniform scalar quantizer considering a uniform
%%    distributed signal
%%    hint: use MSE as distortion metric and plot the SNR
M = [4,8,16,32,64,128];
%cioè faccio le diverse casistiche di num di livelli, cioè con diversi num di bit
%(e quindi di livelli)(length(M)), e per le due distrib

max_s_g = max(s_g);
min_s_g = min(s_g);

max_s_u = max(s_u);
min_s_u = min(s_u);


mse_g = zeros(length(M),1);
mse_u = zeros(length(M),1);

% for each level
for m = 1:length(M)

    % compute the quantization step
    delta_g = (max_s_g-min_s_g) / M(m);
    delta_u = (max_s_u-min_s_u) / M(m);

    % perform quantization
    %%posso usare la mid rise o la mid tread, noi mid tread
    %%procedimento standard
    y_g = round(s_g/delta_g)*delta_g;
    y_u = round(s_u/delta_u)*delta_u;

    y_g_vals = unique(y_g);
    y_g(y_g == y_g_vals(end)) = y_g_vals(end-1);
    y_u_vals = unique(y_u);
    y_u(y_u == y_u_vals(end)) = y_u_vals(end-1);

    % quantization error
    e_g = y_g-s_g;
    e_u = y_u-s_u;

    % MSE
    mse_g(m) = mean(e_g.^2);
    mse_u(m) = mean(e_u.^2);


end

% plot R-D curves
R = log2(M);

snr_g = pow2db(var(s_g) ./ mse_g);
%%pow2db è la stessa cose di 10*log10*
%mse_g è un vettore che contiene gli mse associati ognuno ad una determinato num di bit del quantizers
%quindi anche r_g sarà un vettore
snr_u = pow2db(var(s_u) ./ mse_u);

% Theoretical bound
snr_t = pow2db(2.^(2*R));
%%RIVEDI FORMULA

figure();
plot(R,[snr_g, snr_u, snr_t']);
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');
legend('Gaussian','Uniform','Theoretical Bound');
grid on;

%% 3) Design an optimum uniform quantizer and a Lloyd-Max quantizer for s_g(n)
%%    with the same levels defined in the previous step.
%%    Compare the R-D curves obtained with those quantizers
%%    with those obtained with uniform quantizer

mse_g_uo = zeros(length(M),1);
mse_g_lm = zeros(length(M),1);

% for each level
for m = 1:length(M)

    % compute uniform thresholds
    delta_uo = (max(s_g)-min(s_g))/M(m);
    %%trovo la treshold, il nume delle th sono sempre M+1, ma dato che ho inf e - inf ora avrò M-1 th intermedie
    %%faccio un numero dinamico di soglie
    %%e poi faccio la centratura (-FLOOR)
    %nb. fatti così sono equidistanziati (*DELTA)
    th_uo = [-inf delta_uo * ((1:(M(m)-1)) - floor(M(m)/2)) inf];

    % compute optimum reconstruction levels
    cb_uo = zeros(M(m),1);
    for level = 1:M(m)
        mask = s_g >= th_uo(level) &  s_g < th_uo(level+1);
        cb_uo(level) = mean(s_g(mask));
    end

    % quantize
    s_g_uo = zeros(size(s_g));
    for level = 1:M(m)
        mask = s_g >= th_uo(level) &  s_g < th_uo(level+1);
        s_g_uo(mask) = cb_uo(level);
    end

    % quantization error
    e_g_uo = s_g_uo-s_g;

    % MSE
    mse_g_uo(m) = mean(e_g_uo.^2);

    % Lloyd-Max quantizer with Matlab implementation
    [th_lm,cb_lm]= lloyds(s_g,M(m));
    [~,s_g_lm] = quantiz(s_g,th_lm,cb_lm);

    % quantization error
    e_g_lm = s_g_lm'-s_g;

    % MSE
    mse_g_lm(m) = mean(e_g_lm.^2);


end

% compute SNR
snr_g_uo = pow2db(var(s_g) ./ mse_g_uo);
snr_g_lm = pow2db(var(s_g) ./ mse_g_lm);

% plot R-D curves
figure();
plot(R,[snr_g,snr_g_uo, snr_g_lm, snr_t']);
xlabel('Rate [bit/symbol]');
ylabel('SNR [dB]');
legend('Uniform','Optimum uniform','Lloyd-Max','Theoretical Bound');
grid on;
