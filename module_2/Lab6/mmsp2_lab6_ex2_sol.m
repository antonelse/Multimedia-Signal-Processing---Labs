%% MMSP2 - Lab 6
%  Exercise 2 - ME and MC

clear
close all
clc

%% 1) Load the sequence 'table_tennis.mat' and spatially resize it at half the resolution (hint: use imresize())
load table_tennis.mat

table_tennis = imresize(table_tennis,0.5); %mezza risoluzione significa più schifoso. Twice the resolution è
%l'interpolazione, più sample intermezzi

%% 2) Compute the displaced frame difference (ma intendi SAD??) - use 8x8 blocks, full-search, W=16 pixel, save all motion vectors and visualize them (use quiver())
%rispetto a prima, ora calcolo PER OGNI BLOCCO (e nn solo x uno, come prima) le traslazioni da -16 1 +16 sia in x che in y
[h,w,nf] = size(table_tennis); %nf sono i layer dell terza dimensione, cioè i frame
N = 8;
W = 16;
mv = zeros(h/N,w/N,2);
pred_frame = zeros(h,w);
% loop over all blocks
for r = 1:h/N
    for c = 1:w/N
        % select a block before the next loop !!! Perchè quello che viene dopo (SAD) è riferito ai molteplici confronti di un UNICO blocco
        y0 = (r-1)*N+1;
        x0 = (c-1)*N+1;
        block = table_tennis(y0:y0+N-1,x0:x0+N-1,2);
        % build cost-function
        sad = zeros(2*W+1,2*W+1);
        for dy = -W:W
            for dx = -W:W
                y1 = y0+dy;
                x1 = x0+dx;
                %dato che voglio valutare i motion vector su tutto il possibile range -W:W, dò i limiti della imm generale,
                %dato che se sono agli estremi esco dal boudaries
                %quindi penso ai valori massimi che y1 e x1 possono assumere
                if (y1 < 1 || x1 < 1 || y1+N-1 > h || x1+N-1 > w) %sono i limiti dell'immagine; se esco dall'immagine assegno alla funzione
                  %di costo il valore infinito,e quindi nn prenderò mai
                    sad(dy+W+1,dx+W+1) = +inf;
                else %altrimenti faccio la cosa di prima
                    pred_block = table_tennis(y1:y1+N-1,x1:x1+N-1,1);
                    sad_block = sum(abs(pred_block(:)-block(:))); % SAD ora e' un VALORE, mentre DFD è una matrice di differenza pixel per pixel !!!!
                    sad(dy+W+1,dx+W+1) = sad_block;
                end
            end
        end
        %find minimum, PER OGNI DIVERSA SAD, infatti è dentro i 2 cicli for.
        %OGNI SAD E' RIFERITA AD UN UNICA OPERAZIONE DI TRASLAZ (in x e y) DI UN SINGOLO BLOCCO
        [~,min_idx] = min(sad(:));
        [y1_idx,x1_idx] = ind2sub(size(sad),min_idx);
        y_mv = y1_idx-W-1;
        x_mv = x1_idx-W-1;
        % save motion vectors, che ora me li salvo su una matrice; prima consideravo un solo blocco quindi avevo un solo motion vector di un blocco
        mv(r,c,:) = [y_mv,x_mv]; %per storare le coordinate di movimento x e y come posso fare ? prendo una matrice con 2 layer e sul primo salvo
        %le y e sul secondo layer (il sottostante, pensala il 3d) salvo le x
        % store predicted block into predicted frame
        %r e c sono fissare dentro questi 2 for, quindi lo store avviene nella terza dimensione, potevo fare cioè
        %mv(r,c,1) = y_mv;
        %mv(r,c,2) = x_mv;
        pred_frame(y0:y0+N-1,x0:x0+N-1) = table_tennis(y0+y_mv:y0+y_mv+N-1,x0+x_mv:x0+x_mv+N-1,1);
        %y0 e x0 sono riferiti al blocco che ho scelto(in (:,:,2));
        %applicando questi + i motion vector al reference frame (cioè (:,:,1)), sto dicendo come il blocco in (:,:,2) si evolve  (all'indietro) in (:,:,1);
        %pred_frame lo scandisco lo stesso size del block, infatti con esso ricostruisco tutto il frame, stimandolo
        %attraverso frame 1 predico il frame 2, applicando il motion vector (ricavati prima) ai blocchi del frame 1
    end
end

figure();
subplot(1,3,1);
imagesc(table_tennis(:,:,1),[0 255]);
%[0 255] indicano il min e il max della scala di grigi. Facendo cosi è la stessa
%cosa di fare imshow, vero ??? o dipende calla colorscale gray ?
title('frame 1');
colormap gray;
axis image;

subplot(1,3,2);
imagesc(table_tennis(:,:,2),[0 255]);
title('frame 2');
colormap gray;
axis image;

subplot(1,3,3);
imagesc(pred_frame,[0 255]);
title('frame 2 predicted by mv');
colormap gray;
axis image;

%% compute DFD
dfd = table_tennis(:,:,2) - pred_frame;
%quella di prima era riferita al singolo blocco

% display motion vector superimposed to the frame
figure();
imagesc(dfd,[-255 255]);
colormap gray;
hold on;
x = N/2+1:N:w;
y = N/2+1:N:h; %N/2 poichè li plotto con origine al centro del blocco
quiver(x,y,mv(:,:,2),mv(:,:,1));
%plotta i vettorini

%% 3) Compute mean and variance of DFD and normal frame difference
fd = table_tennis(:,:,2) - table_tennis(:,:,1); %normale differenza tra frame

mean_dfd = mean(dfd(:));
var_dfd = var(dfd(:));

mean_fd = mean(fd(:));
var_fd = var(fd(:));

fprintf('DFD mean: %.2f - var:%.2f\n',mean_dfd,var_dfd);
fprintf('FD mean: %.2f - var:%.2f\n',mean_fd,var_fd);

%% 4) Display DFD and frame difference

figure();
subplot(1,2,1);
imagesc(dfd,[-255 255]);
colormap gray;
axis image;
title('DFD');

subplot(1,2,2);
imagesc(fd,[-255 255]);
colormap gray;
axis image;
title('FD');
