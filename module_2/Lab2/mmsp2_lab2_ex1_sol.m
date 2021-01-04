%% MMSP2 - Lab 2
%  Exercise 1 - Basic scalar Quantization

clc
clearvars
close all

%% Generate 1000 samples with gaussian distribution
rng(21);

x_var = 3;
x = randn(1000,1) * sqrt(x_var);
%rand genera una distr gaussiana, e moltiplicato per la dev stand allargo la campana
%% Quantize with a scalar mid-tread quantizer with fixed quantization step (la x, decision reg) delta=2

delta = 2;

%%%MID TREAD (reproduction region in zero)%%%
y1 = floor(x/delta+0.5)*delta; %equivalent to round(x/delta)*delta;
%%SCALAR MID TREAD: scalar uniforme, NON MIDRISE, cioè nn ha rise in x=0;
%%QUI NN DA I LIVELLI(CIOE' I BIT)

figure();
plot(x,y1,'.');
xlabel('in');
ylabel('out');
grid on;
title('Mid-tread');

%% Quantize with a scalar mid-rise quantizer with fixed quantization step delta=2

delta = 2;

%%%MID RISE (repro contiene lo zero, aka ho il gradino in zero)%%%
y2 = (floor(x/delta) + 0.5)*delta; % equivalent to (round(x/delta-0.5)+0.5)*delta;

figure();
plot(x,y2,'.');
xlabel('in');
ylabel('out');
grid on;
title('Mid-rise');

%% Quantize with a scalar mid-tread quantizer with M=4 output levels
%%qui ho il NUM DI BIT

M = 4;
delta = (max(x)-min(x))/M;

y3 = floor(x/delta + 0.5) * delta; %equivalent to round(x/delta)*delta;

% be aware of the boundaries !
y3_values = unique(y3);
%mappo gli out of bound nell'ultimo valore li mappo al penutimo
%y3_values(end) è il val più grande di y3
%y3 == y3_values(end) sono tutti i val più grandi di y3
y3(y3 == y3_values(end)) = y3_values(end-1);

figure();
plot(x,y3,'.');
xlabel('in');
ylabel('out');
grid on;
title('M = 4');

%% Quantize using cb = [-5,-3,-1,0,1,3,5] as reproduction levels
% and th = [-4,-2,-0.5,0.5,2,4] as thresholds
%cb sono i gradini della y (i livelli)
%th sono i valori di soglia in x

cb = [-5,-3,-1,0,1,3,5];
th = [-4,-2,-0.5,0.5,2,4];

th = [-inf, th, inf];

%costruisco la maschera (cioè regola di quantizzazione) poichè non ho il delta. La creo per ogni intervallo
y4 = zeros(size(x));
for level = 1:length(cb)
    mask = x >= th(level) &  x < th(level+1);
    y4(mask) = cb(level);
    %%il val di y per le x (che sono le mask) corrisponde al livello dello stesso indice
end

figure();
plot(x,y4,'.');
xlabel('in');
ylabel('out');
grid on;
title('Custom');

%% Power and var
e = y1 - x;

Px = mean(x.^2);
Pe = mean(e.^2);
%%è la stessa cosa della varianza per valori che tendono ad infinito

sig2x = var(x);
sig2e = var(e);

snr_p = Px / Pe;
snr_s = sig2x / sig2e;


%% Compute distortion using MSE for each one of the above quantizers
%aka la Pe di sopra

mse1 = mean((y1-x).^2);
mse2 = mean((y2-x).^2);
mse3 = mean((y3-x).^2);
mse4 = mean((y4-x).^2);

%% Compute SNR for each one of the above quantizers

snr1 = var(x) / mse1;
snr2 = var(x) / mse2;
snr3 = var(x) / mse3;
snr4 = var(x) / mse4;

% if we want the snr to be in dB we can use one of the following:
snr1_db = 10*log10(snr1);
snr2_db = db(snr2,'power');
