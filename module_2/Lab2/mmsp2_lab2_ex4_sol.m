%% MMSP2 - Lab 2
%  Exercise 4 - Vector quantization

clearvars
close all
clc

%% 1) Consider a 2D vector quantizer with codebook y1=(1,2), y2=(1,4), y3=(-1,2), y4=(0,-2)
%%    Show optimal assignement regions

%riscrivo il codebook in vettori colonna, quindi prima riga le x e seconda riga le y
cb = [ 1,1,-1,0;
       2,4,2,-2 ];

figure()
%% voronoi, crea il luogo dei punti equidistanti da quei 2 punti (equiprobab)
%creo la griglia voronoi dai punti che ho
voronoi(cb(1,:),cb(2,:)); %(vettore contenenti tutte le x, vettore contenente tutte le y)
hold on;
% plotto i punti che ho
plot(cb(1,:),cb(2,:),'o');
grid on;
%con voronoi, i miei punti saranno per costruzione i centroidi delle regioni
%%semplice "zoom", simil estendi di autocad, sposto min più a sx e max più a dx
xlim([min(cb(1,:))-1,max(cb(1,:))+1]);
ylim([min(cb(2,:))-1,max(cb(2,:))+1]);


%% 2) Quantize the sequence x=(-4:5) using groups of 2 samples at time

x = -4:5;
xgroup = [x(1:2:end) ; x(2:2:end);];
%scrivendo così ho la coppia [x1 ; x2] ottenuta prendendo i sample da x
%e saltando (in modo sfalsato tra i 2) un sample si e uno no

%alloco lo spazio per mettere la x (presa a coppie) quantizzata
xgroup_q = zeros(size(xgroup));

%faccio la quantizzazione
for sample_idx = 1:size(xgroup,2)
    %2 indica che voglio sapere il num di colonne di xgroup, è come lenght ma applicato alla matr
    %scandisco cioè i samples composti dalla coppia x1,x2
    sample = xgroup(:,sample_idx);
    %nb, i sample sono coppie di valori messi a colonna, quindi ognuno sarà indicato da una propria colonna
    %avente indice sample_idx

    % compute the euclidean distance from each centroid, poichè cb sarà il valore che darò al sample,
    %quindi una volta dopo la quantizzazione
    dist = sum((sample - cb) .^2); %ma  SQRT ??? E' l'ERRORE

    % find the nearest centroid
    % min mi restituisce il valore min di dist, ma facendo [~,min_index] richiedo l'INDICE
    % di questo valore minimo del vettore dist (che è in ~, infatti non mi interessa), e lo salvo su min_idx !!!!
    [~,min_idx] = min(dist);
    % scelgo il sample del mio codebook in accordo con l'indice della distanza minima,
    %dato che l'indice di dist ha corrispondenza con quello del codebook, poichè da esso calcolato in linea 49
    sample_code = cb(:,min_idx);
    % fillo la matrice xgroup_q con i valori trovati, quindi quantizzati
    xgroup_q(:,sample_idx) = sample_code;

end

% MATLAB way
vqenc = dsp.VectorQuantizerEncoder('Codebook',cb);
xgroup_q_ml_idxs = vqenc.step(xgroup);
xgroup_q_ml = cb(:,xgroup_q_ml_idxs+1);

%% 3) Plot the original and the quantized sequences

figure();
voronoi(cb(1,:),cb(2,:));
hold on;
leg = [];
leg(1) = plot(xgroup(1,:),xgroup(2,:),'x');
% leg(2) nn funge poichè è della parte del matlab way
leg(2) = plot(xgroup_q_ml(1,:),xgroup_q_ml(2,:),'o');
legend(leg,'Original','Quantized');
grid on;
