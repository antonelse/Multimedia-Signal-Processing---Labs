%% MMSP2 - Lab 6
%  Exercise 1 - ME and MC

clear
close all
clc

%% 1) Load the sequence 'table_tennis.mat' consisting of two grayscale frames

load table_tennis.mat %aka "video" di 2 frame

%% 2) Select the 8x8 block starting at (x,y)=(35,150) of the second frame.
% Perform ME using W=16 pixels and the first frame as reference
N = 8;
W = 16;
y0 = 35;
x0 = 150;
% select the 8x8 pixels block
block = table_tennis(y0:y0+N-1,x0:x0+N-1,2); %prendo una 8x8 del secondo frame (aka secondo layer)
% build cost function testing all possible motion vectors
sad = zeros(2*W+1,2*W+1);

for dy = -W:W
    for dx = -W:W
        pred_block = table_tennis(y0+dy:y0+dy+N-1,x0+dx:x0+dx+N-1,1);
        %prendo il blocchetto 8x8 dal primo frame, spostandomi da -16 px a 16 px
        sad_block = sum(abs(pred_block(:)-block(:)));
        %faccio l'errore di pred tra il frame 1 e il 2, tra le matrici ma unfoldate. Facendo
        %sum diventa un numero, che sarà quindi la mia funzione di costo
        %DIFF TRA FRAME REFERENCE(che traslo) e QUELLO ATTUALE
        sad(dy+W+1,dx+W+1) = sad_block;
        %quindi dy+W+1
        %quindi la matrice sad contiene tutte le differenze che ci sono tra il blocco 8x8 del secondo
        %frame (aka block) (che è il mio freame attuale) e quello preso dal primo frame (pred_block) (il precedente, che è il reference)
        %sponstandoci pixel per pixel (da -16 e +16, compreso lo 0) con un blocco 8x8
    end
end

% find the position of the minimum value of the cost function
[~,idx_min] = min(sad(:));
[y1_idx,x1_idx] = ind2sub(size(sad),idx_min);
%ind2sub trova l'indice del minimo valore della matrice sad !!! modalità difficile ! comoda per trovare indici di più valori in ingr;
%size(sad) serve solo per indicare la struttura della matrice
%è la stessa cosa di fare [index_y,index_x]=min(sad) se c'è UN SOLO min. Se i minimi fossero di più, conviene passargli un vettore
%di minimi e ind2sub ci restituirebbe un vettore per le coordinate x dei minimi e un vettore y con le coordinate y degli stessi minimi
%costruisco ora i motion vector, considerando gli spostamenti
y_mv = y1_idx-W-1; %questi valori indicano di quanto si è spostato il pixel in x e in y dal primo (0) al secondo frame (1)
x_mv = x1_idx-W-1; %tolgo w poichè la matrice SAD era 2W per comodità

%% 3) Display the cost function, the starting block and its estimate
figure();
imagesc(sad);% uses the full range of colors in the colormap
axis image;
colorbar;

%DIFFERENZA TRA IMSHOW E IMAGESC ????? AIDONNNOOOUUUU

figure();
imshow(table_tennis(:,:,1)/255); %uses the default display range for the image data type and optimizes figure
%mostro il primo frame, cioè il reference
%se divido per un num grande, abbasso il valore del pixel, avvicinandomi al nero puro !
hold on;
%ora costruisco il quadrato rosso del quale traccerò lo spostamento. La posiz finale sarà indicata nel secondo(ed attuale) frame con il quadrato verde(vedi sotto)
plot([x0 x0+N-1],[y0 y0],'LineWidth',2,'Color','red');
plot([x0 x0],[y0 y0+N-1],'LineWidth',2,'Color','red');
plot([x0+N-1 x0+N-1],[y0 y0+N-1],'LineWidth',2,'Color','red');
plot([x0 x0+N-1],[y0+N-1 y0+N-1],'LineWidth',2,'Color','red');

x1 = x0+x_mv; %1 è riferito al frame attuale, 0 al precedente (reference frame)
y1 = y0+y_mv;

plot([x1 x1+N-1],[y1 y1],'LineWidth',2,'Color','green');
plot([x1 x1],[y1 y1+N-1],'LineWidth',2,'Color','green');
plot([x1+N-1 x1+N-1],[y1 y1+N-1],'LineWidth',2,'Color','green');
plot([x1 x1+N-1],[y1+N-1 y1+N-1],'LineWidth',2,'Color','green');
