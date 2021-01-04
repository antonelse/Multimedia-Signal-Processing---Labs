clc, clearvars, close all


%% 1) Load the file 'voiced_a.wav' and consider only a 300ms frame.
% Plot the magnitude of the frequency response of the frame

[s,Fs] = audioread('voiced_a.wav');

dur = 0.3;
s = s(1:dur*Fs,:); %entrambi i channel L e R, è una matrice

% s = sin(2*pi*((1:dur*Fs)-1)*101/Fs);  % just a check that a pure tone exhibits an ordered set of peaks at L, 2L, 3L, etc.
%%BHO ???

N = size(s,1);
S = fft(s);

% faxes = (0:N/2-1)/N*Fs; % asse x delle freq, fft centrate in zero, bilatero
faxes = (0:N-1)/N*Fs; %fft tutta a dx (da 0 a N-1, con intermezzo 1/Fs); *Fs è al num quindi SONO IN FREQ;
%al denominatore  sarebbe stato un asse temporale

spectrum_fig = figure();
%semilogx(faxes,db(abs(S(1:floor(end/2)))),'DisplayName','S');
semilogx(faxes,db(abs(S)),'DisplayName','S'); %faccio l'asse logaritmico SOLO in x, cioè asse delle freq,
% come scritto sul comando


xlabel('f [Hz]');
ylabel('Magnitude [dB]');
% xlim([0,Fs/2]);
xlim([0,Fs]); %setta limite max e minimo nella stampa dell'asse x. Prendo a 0 a Fs
grid on;
legend();

% pause;
%% 2) Perform pitch detection using auto-correlation method.
%  Consider only frequencies between 60 Hz and 500 Hz


% compute autocorrelation (hint: use two output arguments of xcorr() and
% pay attention to MATLAB normalization)
[r,rlags] = xcorr(s,s,'coeff'); %faccio l'autocorr e rlags è la Tao; r NON e' la matrice di corr,
%ma SOLO un vettore che contiene i coefficienti della stessa, quindi rpos e rc sono anch'essi vettori
%'normalized' or 'coeff' — Normalizes the sequence so that the autocorrelations at zero lag equal 1

% consider correlation only for positive lags, including 0
rpos = r(rlags >=0); %segnale reale, quindi autocorr simmetrica
%qui prendo i valori di autocorr associati a lag positivi (che essendo reali sono simmetrici con i lag negativi)

r_fig = figure();
plot(rlags(rlags>=0),rpos); %plotto l'autocorr positiva, calcolandomi asse x e asse y
title('r');
xlabel('lag');
grid on;

% find maximum peak within the accepted range
Fmin = 60;
Fmax = 500;
lagmin = floor(Fs/Fmax); %è un numero puro
lagmax = ceil(Fs/Fmin);

rc = rpos(lagmin+1:lagmax+1);
%+1 needed to index coherently with matlab
%prendo così i coeff associati ai lag d'interesse, cioè quelli del range di freq da trovare

[rc_maxv,rc_maxi] = max(rc); %mi ritorna l'indice(rc_maxi) del valore massimo del
%vettore rc nonchè il valore stesso(v sta per value)
%rc_maxi is indexed in matlab way (starting from 1)

r_maxi = lagmin+rc_maxi-1; %definisco questi vettore perchè:
%ho trovato il massimo dell'autocorr restringendo il campo tra le 2 freq (vedi rc),
%poi però necessito del lag totale quindi devo riconsiderare tutti  coeff dell'autocorrelazione(r)
%riaggiungendo l'offset iniziale (IN INDICE e non in valore) imposto a rc (cioè a rpos) il quale è lagmin
%rc_maxi-1 compensates for Matlab indexing, so that rc_maxi is 0 based indexed

r_maxv = rc_maxv; %il valore è lo stesso. Avevo solo problemi di INDICE

% hold on;
% stem(r_maxi,r_maxv,'r');

pitch_lag = r_maxi; %quindi ho l'indice temporale (in samples) associato al termine di autocorr più alta,
%cioè al valore freq più definito
pitch = Fs/pitch_lag; % è la mia freq ottenuta come ripetizione del segnale (autocorr method)

fprintf('Pitch lag: %d samples - freq: %.2f Hz\n',pitch_lag,pitch); %stampa a linea di comando

% pause;
%% 3) Compute LPC coefficients of order 12
% hint: if using correlation, consider only positive lag values.
% Use toeplitz() to build correlation matrix from correlation vector
%quindi mi costruisco la matr di corr partendo dai coeff contenuti in r

p = 12; %l'ordine è il lag max del sist
R = toeplitz(rpos(1:p)); %coefficients from 0 to p-1
%prende il vettore r (solo la parte positiva, come scritto nell'hint, e di lungh uguale all'ordine) e
%lo pone come prima riga della toeplitz R. Da qui costruisce le diag costanti (easy)
%toeplitz: simm a diag costanti
a = R\rpos(2:p+1); %rpos(2:p+1) : coefficients from 1 to p
%BACKSLASH; risolvo il sistema R*a = rpos(), dato che mi devo trovare coeff a del sist LPC
%NON prendo il valore zero poichè dalla def di lpc ( vedi pag 13 del pdf esercitaz )

% Alternatively, you can use the lpc() function, be aware of what is returned
%a_lpc = lpc(s,p);

lpc_fig = figure();
stem(a);
title('a');
xlabel('idx');
ylabel('Value');
grid on;

%in questo passaggio calcolo la traformata SOLAMENTE del DENOMINATORE (compresa la parte 1-, aggiunta in testa) e non della
%fdt H(z) ( vedi pag 13 del pdf esercitaz ). La fdt totale la scrive direttamente sul plot sotto INFAME !
A = fft([1 ;-a],N);

figure(spectrum_fig); %riga 21 !! aggiunge altro grafico alla figura spectrum_fig definita nella RIGA 21
hold on;
% plot(faxes,db(abs(1./A(1:floor(end/2)))),'DisplayName','1/A');
plot(faxes,db(abs(1./A)),'DisplayName','1/A');
%(1./A) è lo SHAPING FILTER !!!!!!!!!

legend();

% pause;
%% 4) Plot the prediction error and its magnitude spectrum

e = filter([1; -a],1,s);
%alimentando l'inverso dello shaping (quindi solo A(z)) con il segnale s, ottengo l'errore !!!
%vedi il diagramma (o la trasf in zeta) dell LPC ( vedi pag 13 del pdf esercitaz )

time_fig = figure();
plot((0:length(e)-1)/Fs,e,'DisplayName','e'); % /Fs quindi è asse temporale
xlabel('Time [s]');
ylabel('Value');
grid on;

E = fft(e);
figure(spectrum_fig);
hold on;
% plot(faxes,db(abs(E(1:floor(end/2)))),'DisplayName','E');
plot(faxes,db(abs(E)),'DisplayName','E');

legend();

%% 5) Build an impulse train with the estimated pitch period
% hint: initialize a vector of zeros and put a 1 every 1/pitch (aka un pitch_lag) seconds

u = zeros(length(s),1);
u(1:pitch_lag:end) = 1; % è l'hint, questo lag corrisponde al tempo della freq pitch (1/pitch è un tempo)

% normalize the energy: force the energy of the
% impulse train to be equal to that of the residual
u = u / std(u) * std(e); % energy normalization
%prima lo normalizzo a 1 e poi alla dinamica di e

figure(time_fig); %le plotto sempre sul grafico time_fig, aka riga 125
hold on;
plot((0:length(u)-1)/Fs,u,'DisplayName','u');
legend();

% pause;
%% 6) Consider the impulse train as excitation and build synthetic speech
%hint: if curious use fvtool() to visualize the filter you are using.
% Can you see the formants?

s_synth = filter(1,[1; -a],u); %faccio il sist a blocchi di pag 13 esercitaz
s_synth = s_synth/max(abs(s_synth)); %normalizzazione del segnale sitetizzato

S_synth = fft(s_synth,N);
figure(spectrum_fig);
hold on;
% plot(faxes,db(abs(S_synth(1:floor(end/2)))),'DisplayName','S synth');
plot(faxes,db(abs(S_synth)),'DisplayName','S synth');

legend();

pause;
%% 7) Listen to the original and the synthetic speech

sound(s,Fs);
pause(1);
sound(s_synth,Fs);
