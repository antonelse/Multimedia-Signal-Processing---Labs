%% MMSP2 - Lab 4
%  Exercise 3 - Vocoder with voiced/unvoiced classification

clear
close all
clc

%% 1) Load the files 'a.wav' and 'shh.wav' and build a single signal
% concatenating them

[s_a,Fs] = audioread('a.wav');
s_shh = audioread('shh.wav');

s = [s_a;s_shh]; %i samples li metto consecutivi in un'unica colonna, entrambi sono vect colonna !!!

%% 2) Loop over every frame and compute voicing detection parameters
% (i.e., zcr and ste) for each one of them

% generate Hamming window
% hint: use hamming()
frame_dur = 0.04; %durata in sec
frame_stride = 0.01; %è il passo delle finestre, aka "overlapp", il release che sarà uguale all'attack (in sec)

frame_len = round(frame_dur * Fs); %durata in samples
frame_step = round(frame_stride*Fs); %durata fade in (cioè fade out)

win = hamming(frame_len);

N = floor((length(s)-frame_len)/frame_step)+1; %lenght(s)/frame_step di indica
%il num di passi che "entrano" nel segnale s cioè il num delle finestre che userò.
%Togliere una lunghezza della finestra mi permette di non "sforare" con l'ultima
%finestra cosi da finire "a pelo" con ultima finestra e segnale
%+1 nn sono sicuro, forse per il floor

% container for classification parameters
parameter = zeros(N, 2);

%analizzo ogni chunks finestrato (analizzo con zcr e ste)
for n=1:N % for each frame, cioè per ogni "finestra"
    %faccio il finestramento, cioè moltiplico chunks di segnale per la win di hamming
    % select a frame and multiply it with the window
    frame = s((n-1)*frame_step+1:(n-1)*frame_step+frame_len).*win;
    %dentro s metto l'indice di s che prenderò. Seleziono il chunks di s che verrà moltiplicato per la win

    % Zero-crossing rate, conto i passaggi per l'asse x, aka zeri
    frame_zcr = sum(abs(diff(frame>0)))/frame_len; % FEDE !!!!!!!! INVAME
    %calculates differences between adjacent elements of X along the first array
    %Y = diff(X) ritorna Y = [X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)]

    % Short-time energy
    frame_ste = sum(frame.^2);

    parameter(n,:) = [frame_zcr,frame_ste];

end

%% 3) Voiced / Unvoiced classification
% Compute thresholds using median()

th = median(parameter);
% https://it.mathworks.com/help/matlab/ref/median.html?searchHighlight=median&s_tid=doc_srchtitle
%nb li devi ORDINARE SORTTTT
%ho preso la mediana dello zcr e ste, quindi è un vect riga di 2 param

% Decision
voiced = parameter(:,1) < th(1) & parameter(:,2) > th(2);
% voiced(i) = 1 if the i-th frame is voiced, 0 if unvoiced
%parameter(:,1) prendo tutto gli zero crossing
%th(1) è la mediana dello zcr
%parameter(:,2) prendo tutto gli ste
%th(2) è la mediana dello ste
%questa è la formula della slide 9 dell'esercitazione;
%per essere voiced si devono verificare entrambe le condizioni


% plot the parameters with the thresholds
figure();
subplot(3,1,1);
plot(parameter(:,1),'r');
hold on;
plot(th(1)*ones(N,1),'r--');
% serve per dare riferimento di tresholf fissa, è la tratteggiata
title('ZCR');
grid on;

subplot(3,1,2);
plot(parameter(:,2),'b');
hold on;
plot(th(2)*ones(N,1),'b--');
title('STE');
grid on;

subplot(3,1,3);
plot(voiced);
title('Voiced');
grid on;
%quindi e voiced quando sono sotto soglia per gli zcr e sopra soglia per ste

%% 4) LPC analysis and synthesi, nb.sto facendo l'inverso, aka feeddo con il segnale e trovo la e
p = 12; % order of the predictor

Fmin = 60;
Fmax = 500;

lagmin = floor(Fs/Fmax);
lagmax = floor(Fs/Fmin);

% container for the synthesized speech
s_synth = zeros(length(s),1);
% s_synth = [];
for n = 1:N

    % select a windowed frame; lo faccio in 2 passaggi, prima indice e poi frame
    frame_idxs = (n-1)*frame_step+1:(n-1)*frame_step+frame_len;
    frame = s(frame_idxs).*win;


    %% 4a) Compute LP coefficients and prediction error

    [r,rlags] = xcorr(frame,frame,lagmax); % limits the lag range from -maxlag to maxlag for either of the previous syntaxes
    rpos = r(rlags >=0);

    R = toeplitz(rpos(1:p));
    a = R\rpos(2:p+1); %risoluz sist lin Ax=B

    if n == 1 %faccio l'if per motivi di predizione del primo frame
        [e, Sf_e] = filter([1; - a], 1, frame); %Sf_e è la z-trasf del filtro, sto facendo la STFT
    else
        [e, Sf_e] = filter([1;- a], 1, frame, Sf_e); % Sf_e agisce da initial conditions, è la trasf poichè necessito
        %di una condi già trasf
    end


    %% 4b-1) Voiced segment: %gli n sono i frames
    if voiced(n) == 1
        % Pitch detection
        rc = rpos(lagmin:lagmax); %focalizzandoci su un range limitato di freq, imposto i lag relativi a queste

        [~,rc_maxi] = max(rc); %rc_maxi è l'indice del massimo, nn mi interessa il valore
        pitch_lag = lagmin+rc_maxi; %rimuovo l'offset

        % Generate impulse train
        u = zeros(length(frame),1);
        u(1:pitch_lag:end) = 1;


    %% 4b-2) Unvoiced segment:
    else
        % Generate random noise instead of impulse train
        %hint: use randn()
        u = randn(length(frame),1); %randn ritorna una matrice, come zeros, ma impostando a
        %1 la seconda dim,ho un vett colonna

    end

    %% 4c) Normalize the energy of the excitation signal
    % (i.e., make it equal to that of the prediction error)
    u = u / std(u) * std(e);

    %% 5) Shaping filter (i.e., synthesize the frame) %ora mando la u a feddare il filtro quindi sinstetizzo s
    %hint: pay attention to the overlap
    if n == 1 %come quello di prima ma schema visto al contrario, fase di sintesi
        [frame_synth, Sf_x] = filter(1, [1; -a], u);
    else
        [frame_synth, Sf_x] = filter(1, [1; -a], u, Sf_x);
    end

    s_synth(frame_idxs) = s_synth(frame_idxs) + frame_synth; %faccio l'OLA, e quindi sommo il frame_synth calcolato
    %all'iterazione n con il "frame" s_synth dell'iterazione precedente ma shiftato dello step !!! quindi OLA

end


%% 6) Listen and compare the original and synthetic signals

% normalize in [-1;1] range just for listening
s_synth = (s_synth - min(s_synth))/(max(s_synth) - min(s_synth))*2 -1; %range +1 e -1

% sound(s,Fs);
% pause(5);
% sound(s_synth,Fs);
