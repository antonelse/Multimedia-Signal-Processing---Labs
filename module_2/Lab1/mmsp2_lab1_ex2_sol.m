%% MMSP2 - Lab 1
%  Exercise 2 - Image signal encoding

clearvars
close all
clc

%% 1) Load the image 'lena512color.tiff' and display the normalized histogram
%%    for all (R,G,B) components
%%    hint: 'lena512color.tiff' is an 8-bit RGB image, better to convert into double
%%    hint: express the RGB components as vectors

im = imread('lena512color.tiff');
%imread restituisce una matrice a 3 dimensioni con righe-colonne
%(2D part) data dalla dimensione dell'immagine in pixel
%e 3 layer (la terza dimensione) che rappresentano rispettivamente Rosso Green Blue
disp(['size(im): ' mat2str(size(im))]);

R = double(im(:,:,1));
G = double(im(:,:,2));
B = double(im(:,:,3));

%prendo TUTTE le righe di R e le metto una di seguito all'altra
%(512x512 è la LUNGHEZZA del singolo vettore colonna R poichè dispongo 512 righe
%lunghe ognuna 512 elementi una di seguito all'altra)
%questo mega vettore colonna contiene i valori di Rosso di ogni pixel dell'immagine
R = R(:);
G = G(:);
B = B(:);

alphabet = 0:255;
d_R = hist(R, alphabet);
d_G = hist(G, alphabet);
d_B = hist(B, alphabet);
p_R = d_R/sum(d_R);
p_G = d_G/sum(d_G);
p_B = d_B/sum(d_B);

figure()
colors = {'red','green','blue'};
bar(alphabet,p_R,'red','FaceAlpha',0.7);
%aka(numero barre,altezza barre,colore barre,trasparenza a 0.7)
hold on; %per permettere di plottarle tutte in un unico grafico, figura
bar(alphabet,p_G,'green','FaceAlpha',0.7);
bar(alphabet,p_B,'blue','FaceAlpha',0.7);
grid('on');

%% 2) Compute the entropy of each channel

H_R = -sum(p_R(d_R > 0) .* log2(p_R(d_R > 0)));
H_G = -sum(p_G(d_G > 0) .* log2(p_G(d_G > 0)));
H_B = -sum(p_B(d_B > 0) .* log2(p_B(d_B > 0)));

disp(['H Red channel = ', num2str(H_R)]);
disp(['H Green channel = ', num2str(H_G)]);
disp(['H Blue channel = ', num2str(H_B)]);

%% 3) Let X represent the source of the red channel and Y the source of the
%%    green channel. Compute and show p(X,Y).
%%    hint: use the function imagesc() to show p(X,Y)

d_joint = hist3([R G], {alphabet, alphabet});
disp(['size(d_joint): ' mat2str(size(d_joint))]);
%Nb: è una matrice
p_joint = d_joint/sum(d_joint(:));
%Nb: è ancora una matrice, avendo fatto una semplice divisione per uno scalare

figure();
imagesc(db(p_joint)), axis xy;
%indico a quali assi corrispondono i 2 valori ricavati da imagesc
title('p(X,Y)');
ylabel('X symbols');
xlabel('Y symbols');

%% 4) Compute and display the joint entropy H(X,Y) and verify that H(X,Y) ≤ H(X) + H(Y).
% NON CONFONDERE JOINT ENTROPY CON CONDITIONAL ENTROPY

H_joint = -sum(sum(p_joint(d_joint > 0) .* log2(p_joint(d_joint > 0))));
fprintf('joint entropy: %.3f bit/2 pixel <= %.3f bit/2 pixel\n',H_joint, H_R + H_G);


%% 5) Suppose to encode Y (aka green channel) with H(Y) bits and transmit source N = aX+b-Y instead of X
%%    Compute the entropy of N and the conditional entropy H(X|Y)
% Compute coefficients a and b with LS (LEAST SQUARES !!!!)

X = R;
Y = G;

%per trovare i coeff impongo la eq a zero, cioè aX+b-Y=0
coeff = [X ones(size(X,1),1)]\Y;
%(si chiama mldivide) è x = A\B  and solves the system of linear equations A*x = B (A*coeff=Y)
%con A=[X ones(size(X,1),1)] nel nostro caso, dovendo considerare i 2 coeff
%The matrices A and B must have the same number of rows.
%X e Y sono entrambi vet colonna della stessa dim
%[X ones(size(X,1),1)] significa fare una matrice (nota le quadre) con prima colonna X
%e seconda colonna (solo una, poichè il numero di colonne è specificato dal secondo elemento di ones) composta da
%tutti "1", lunga quanto X (vedi il size(x,1))
%la seconda colonna tiene conto del coeff di b, che è 1 quindi nella colonna ci saranno tutti "1"
%coeff è un vett colonna a 2 righe (è una 2x1) poichè facendo la divisione per Y torna così !!!
% A è una lenght(X)*2 e coeff è una 2*1 quindi si può fare
disp(['size(coeff): ' mat2str(size(coeff))]);

a = coeff(1);
b = coeff(2);

N = round(a*X+b-Y);
%nb: N è ancora un vettore colonna della stessa dim di X nonche di Y

% Compute H(N)
alphabet = unique(N);
%elimina le occorrenze uguali (no ripetiz) di un vettore e le ordina in senso crescente (da 0 a 255)
%genero l'alfabeto partendo dai valori che ho, cioè quelli di N
d_N = hist(N(:),alphabet);
p_N = d_N/sum(d_N);
H_N = -sum(p_N(d_N > 0) .* log2(p_N(d_N > 0)));
fprintf('entropy N: %.3f bit/pixel\n',H_N);

% Conditional entropy H(X|Y) = H(X,Y) - H(Y)
H_cond_X = H_joint - H_G;
%G poichè il green è la Y
fprintf('cond entropy X|Y: %.3f bit/pixel\n',H_cond_X);
