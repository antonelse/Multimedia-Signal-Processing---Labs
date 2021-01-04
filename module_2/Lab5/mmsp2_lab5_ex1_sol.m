%% MMSP2 - Lab 5
%  Exercise 1 - DCT/KLT comparison

clear
close all
clc


%% 1) Load the image 'lena512color.tiff' and extract the luminance component
im = imread('lena512color.tiff');
R = double(im(:,:,1));
G = double(im(:,:,2));
B = double(im(:,:,3));

%luminance è CL dell'RGB attaverso una matrice di trasformazione
Y = 0.299*R + 0.587*G + 0.114*B;

[Nr,Nc] = size(Y); %scritto così, restituisce num di righe Nr e num di colonne Nc
%% 2) Consider the first 8x8 pixels block and compute its 8x8 DCT coefficients.
% Use different methods and compare them.

N = 8; % dimension of the image block
K = 8; % dimension of the projection space, k-esime freq, basi della DCT

block = Y(1:N,1:N); %prendo la luminace dell'imm, solo 8x8
block_dct = zeros(K,K,4);
% 4 è la lunghezza della terza dimensione, mi servirà per il confronto dei metodi di calcolo di dct

% Method 1: use the separability property of DCT
% Define transform matrix using the given equation
% Apply transform matrix by rows and by columns


Tdct1 = zeros(K,K);
%scrivo la formula che ho sulle slide del lab
for k = 1:K %baso k sul num delle freq discrete che descriveranno la mia imm
    if k == 1
        a = sqrt(1/N);
    end
    if k ~= 1 %diverso ~=
        a = sqrt(2/N);
    end
%prima ho definito il coeff che moltiplica e poi faccio con il cos
    for l = 1:N %baso l sul num dei pixel %MA COSA INDICA DI PRECISO ???
        Tdct1(k,l) = a*cos((2*l-1)*pi*(k-1)/(2*N));
    end

end

%CALCOLO I COEFF, che sono contenuti in block_dct
block_dct(:,:,1) = Tdct1 * block * Tdct1'; %prima dimensione 3D (aka primo layer).
%Faccio la formula delle slide sulla separabilità

% Method 2: as method 1 but using the function dctmtx() to build the transform matrix
Tdct2 = dctmtx(K); %costruisco la trasf con K=freq basi
block_dct(:,:,2) = Tdct2 * block * Tdct2'; %seconda dimensione 3D (aka secondo layer).
%stessa cosa di sopra

% Method 3: use function dct2()
block_dct(:,:,3) = dct2(block); %terza dimensione 3D (aka terzo layer).
%SONO GIA' I COEFF. Prima era la matrice di trasf. QUI FACCIO DIRETTAMENTE LA TRASFORMAZIONE

% Method 4: define the transformation g(n1,n2,k1,k2) and apply it, cioè facciamo la DCT in 2D di sopra ma a mano
% vedi slide lab x le formule
g = zeros(N,N,K,K);
for n1 = 1:N
    for n2 = 1:N
        for k1=1:K
            for k2=1:K
                a_k1 = sqrt(2/N);
                if k1 == 1 %sarebbe la casisitica di k=0, ma l'indicizzazione di matlab fa cagare
                    a_k1 = sqrt(1/N);
                end
                a_k2 = sqrt(2/N);
                if k2 == 1
                    a_k2 = sqrt(1/N);
                end
                g(n1,n2,k1,k2) = a_k1*cos((2*(n1-1)+1)*pi*(k1-1)/(2*N))*a_k2*cos((2*(n2-1)+1)*pi*(k2-1)/(2*N));
                block_dct(k1,k2,4) = block_dct(k1,k2,4) + g(n1,n2,k1,k2)*block(n1,n2);
                %la somma del block_dct iniziale serve per sommare i contributi ottenuti con le iterazioni di k1,
                %poi k2, e poi per ogni n2 e n1, altrimenti senza di esso il fill della matrice viene annullato al ciclo
                %for successivo (a salire) e non considero più i contributi degli annidati; si vede bene eseguendo per
                %step la funzione dentro a questo for
            end
        end
    end
end

%HO FINITO I 4 METODI (4 layers), ora devo quantizzare questi coeff

% Compute the MSE between coefficients in the different cases, quindi in sede di esame nn mi serve.
% all'esame l'mse è dato da
% mean((block(:,:) - block_dct(:,:))^2)
% ed è relativo al segnale vero e proprio

block_mse = zeros(4,4); %devo confrontare ogni layer con tutti gli altri quindi ho 4x4=16 combinazioni con le quali fillerò
%la mia matrice che mi indicherà i valori mse di ogni singola combinazione
for idx1 = 1:4
    for idx2 = 1:4
        block_mse(idx1,idx2) = (sum(sum((block_dct(:,:,idx1)-block_dct(:,:,idx2)).^2)))/(N.^2);
        %otterrò una matrice simmetrica; perchè il confronto tra "a" e "b" con quello "b" e "a" è lo stesso ovviamente.
        %Il confronto con se stesso (diagonale) è nullo

    end
end

disp('MSE');
disp(block_mse);

% Are really all these methods equal?
% Yes, the small difference is only due to numerical limitations, cioè sono exp moolto negativi,
%quindi num piccoli, approssimaz matlab

%% 3) JPEG baseline coding - estimate the PSNR (Peak signal-to-noise ratio) and display the reconstructed image
load('Qjpeg.mat'); % quantization matrix
Q = 1; % scaling factor (fine(fine) --> coarse(grossolano)); + è alto e più agisce la quantizz, quindi più fa schifo
Qmatrix = QJPEG * Q; %QJPEG è stata definita automaticamente nella fase di load, è una matrice 8x8

% store the quantized symbols to compute entropy later on of them, so we stored it in layers
[h,w] = size(Y); %[righe(altezza, height),colonne(larghezza, width)]
num_block = h/N*w/N; %quanti blocchi 8x8 entrano nella mia imm per quantizzarla? num_block; sono le DECISION INTERVAL
jpeg_symbols = zeros(K,K,num_block); %ogni layer è un simbolo jpeg, e sono tanti quanti i blocchi; sono le REPRODUCTION REGION
block_idx = 1;

Y_jpeg = zeros(size(Y));
% loop over each 8x8 pixels block (no overlap)
for r = 1:size(Y,1)/8 %size(Y,1)/8 sono il num di blocchi da 8 che entrano nella Y (che è 512x512), in questo caso in una riga (vedi ( ,1));
                      %mi entrano 64 blocchi (in riga, e 64 in colonna) da 8 nella mia imm
    for c = 1:size(Y,2)/8 %r righe, c colonne; in questo caso in una colonna (vedi ( ,2));
        block_r_idxs = (r-1)*N+1:r*N; %parto da 1 come indice(+1) ma considero il blocco di indice 0 (r-1)
        block_c_idxs = (c-1)*N+1:c*N;

        %ora analizzo i 4 step di trasformazione dell'immagine (a,b,c,d)
        % 3a) consider block of size 8x8, prendo una 8x8 all'interno di Y; block_r_idxs è un intervallo di numeri quindi un insieme di righe
        block = Y(block_r_idxs,block_c_idxs);
        % 3b) compute the DCT (function dct2()), quantizzo l'IMMAGINE, e non i coeff. Quelli
        block_dct = dct2(block-127); %voglio valori da -2^(B-1) a 2^(B-1), quindi devo centrarli in zero, abbassando tutti di 2^(B-1)=127
        %l'offset è dettato dallo standard
        %qui con dct2 ho già compresa la trasformazione
        % 3c) threshold quantization, pag 136; la statistica è stata fatta offline ed è contenuta della Qmatrix
        block_dct_coeff = round(block_dct./Qmatrix); %è la quantizzazione dei coeff dell'immagine (o meglio di un suo blocco)
        block_dct_q = block_dct_coeff.*Qmatrix;

        jpeg_symbols(:,:,block_idx) = block_dct_coeff; %sto fillando la terza dim, che mi rappresenta il simbolo(8x8) (cioè i coeff quantizzati)
        %che dò al rispettivo blocco dell'imm
        %OGNI LAYER HA I COEFF ASSOCIATI AD UN BLOCCO. UN LAYER è UN SIMBOLO !!
        %i simboli ci servono solo per calcolarci dopo l'entropia, ma in realtà qui nn mi serve
        block_idx = block_idx + 1;

        % 3d) reconstruct the image from quantized coefficients (function idct2()); sono a questo punto al decoder
        block_idct = idct2(block_dc_q);

        Y_jpeg(block_r_idxs,block_c_idxs) = block_idct + 127; %riaggiungo la centratura post trasformazione
    end
end

% display image and compute PSNR and entropy

% figure()
% subplot(1,2,1);
% imshow(uint8(Y)); %ogni elem della Y è rappresentato con scala di bit a 8-bit unsigned integer.
% title('Original');

% Displays the grayscale image in a figure
% subplot(1,2,2);
% imshow(uint8(Y_jpeg));
% title('JPEG');

mse_jpeg = mean((Y(:)-Y_jpeg(:)).^2);
PSNR_jpeg = pow2db(255.^2/mse_jpeg); %nell'RGB si ha range da 0 a 255, e 255^2 indica la varianza
%massima possibile, che ci serve per il Peak to Noise Ratio

disp(['PSNR JPEG: ' num2str(PSNR_jpeg)])

jpeg_symbols = jpeg_symbols(:); %unfoldo la jpeg_symbols, quindi 64x64(numero dei blocchi, cioè num di simboli)x8x8(valori per ogni simbolo)
jpeg_symbols_values = unique(jpeg_symbols); %alfabeto in ordine, è la roba che facevamo all'inzio del course
jpeg_symbols_count = hist(jpeg_symbols, jpeg_symbols_values);
jpeg_symbols_prob = jpeg_symbols_count./sum(jpeg_symbols_count); %sum(jpeg_symbols_count) sono quindi ora le OCCORRENZE del simbolo
%(cioè è la stessa cosa di size(jpeg_symbols)
jpeg_entropy = -sum(jpeg_symbols_prob.*log2(jpeg_symbols_prob));

fprintf('Entropy JPEG: %.4f bit/symbol\n',jpeg_entropy);

%% 4a) Reconstruct the image using only the DC component (f=0) of the DCT - estimate the PSNR and display the reconstructed image

% store the quantized symbols to compute entropy later on
dc_symbols = zeros(K,K,num_block);
block_idx = 1;

dc_only_mask = false(K,K); %stessa cosa di fare zeros(K,K)
dc_only_mask(1,1) = true; %imposto a 1 solo la freq zero
%faccio cioè per creare una matrice che moltiplicata mi farà agire solo il coeff della dc

Y_dc = zeros(size(Y));
%ora farò uguale a prima (linea 117) ma questo caso le freq basi sono solo la continua
for r = 1:size(Y,1)/8
    for c = 1:size(Y,2)/8

        block_r_idxs = (r-1)*N+1:r*N;
        block_c_idxs = (c-1)*N+1:c*N;

        % 4a) consider block of size 8x8
        block = Y(block_r_idxs,block_c_idxs);

        % 4b) compute the DCT
        block_dct = dct2(block-127);

        % 4c) keep only DC
        block_dc_coeff = round(block_dct./Qmatrix);
        block_dc_coeff(~dc_only_mask) = 0; %CONDIZIONE NOT: ~dc_only_mask; cioè setto a zero
        %gli elementi di block_dc_coeff con indici falsi in dc_only_mask (falsi, cioè = 0);

        %non era la stessa cosa di moltiplicare i coeff per la MASCHERA ???
        %cioè block_dc_coeff = block_dc_coeff.*dc_only_mask

        block_dc_q = block_dc_coeff .* Qmatrix;

        dc_symbols(:,:,block_idx) = block_dc_coeff;
        block_idx = block_idx + 1;

        % 4d) reconstruct the image from quantized coefficients

        block_idct = idct2(block_dc_q);
        Y_dc(block_r_idxs,block_c_idxs) = block_idct + 127;
    end
end

% display image and compute PSNR

figure()
subplot(1,2,1);
imshow(uint8(Y));
title('Original');

subplot(1,2,2);
imshow(uint8(Y_dc));
title('DC');

mse_dc = mean((Y(:)-Y_dc(:)).^2);
PSNR_dc = pow2db(255.^2/mse_dc);

disp(['PSNR DC: ' num2str(PSNR_dc)])

dc_symbols = dc_symbols(:);
dc_symbols_values = unique(dc_symbols);
dc_symbols_count = hist(dc_symbols, dc_symbols_values);
dc_symbols_prob = dc_symbols_count./sum(dc_symbols_count);
dc_entropy = -sum(dc_symbols_prob.*log2(dc_symbols_prob));

fprintf('Entropy DC: %.4f bit/symbol\n',dc_entropy);

% any comment? qui la entropia è più bassa poichè usando solo un coeff,
%èètrasmetto meno bit per simbolo. Ricordo che ogni simbolo è una matr 8x8 associata ad un blocco

%% 4b) Reconstruct the image using only one AC component of the DCT
% Do not quantize, just for the sake of reconstruction

Y_ac = zeros(size(Y));

% fix one component (prendo una "freq")
k1 = 2;
k2 = 3;

for r = 1:size(Y,1)/N
    for c = 1:size(Y,2)/N
        %% 4a) consider block of size 8x8
        block = Y((r-1)*N+1 : r*N, (c-1)*N+1 : c*N);
        %% 4b) compute the DCT
        block_dct = dct2(block);
        %NON FACCIO IL CICLO come prima su k1 e k2 freq poichè ne prendo una già fissata sopra
        %% 4c) keep only coeff (k1, k2)
        coeff = block_dct(k1,k2); %invece di fare una maschera (CHE AVREBBE FUNZIONATO, infatti vedi dopo che fa), qui prendo solo il coeff associato alla mia freq
        block_dct = zeros(size(block_dct)); %rimetto a zero la trasf del blocco, poichè voglio sapere solo il contrib della mia freq
        block_dct(k1,k2) = coeff; %la refillo con il coeffiente di mio interesse, SOLO CON LUI !!!!
        %% 4d) reconstruct the image from quantized coefficients
        %ricostruisco per blocchi, indicando gli indici (prima era la stessa cosa) per Y
        Y_ac((r-1)*N+1 : r*N, (c-1)*N+1 : c*N) = idct2(block_dct);
    end
end

figure()
subplot(1,2,1);
imshow(uint8(Y));
title('Original');

subplot(1,2,2);
imshow(Y_ac, []);
title('AC (2,3)');

%% 5) Consider blocks of dimension 8x8 and estimate the correlation matrix

% compute image blocks
imblocks = im2col(Y,[N,N],'distinct'); %prende dei blocchi NxN da Y, e poi li srotola ponendoli
%in colonna, e ad ogni colonna corrisponde un blocco. Il num di righe è il num degli elementi di ogni blocco.
%Il num delle colonne è il num dei blocchi
%ATTENTO, se nn scrivo 'distinct' fa 'sliding' cioè prende i blocchi overlappati !! quindi avrò un numero altissimo di blocchi (cioè di colonne)

%è la stessa cosa di reshape ma per matrici, invece che per vect !

% compute and remove the mean block
block_mean = mean(imblocks,1);
% mean(A,2) is a COLUMN vector containing the mean of each ROWS
% mean(A,1) is a ROW vector containing the mean of each COLUMNS

imblocks_zm = imblocks - block_mean; %zm aka zero mean, è una matrice contenente i blocchi (in ogni colonna) ma senza val medio

%AUTOCORRELAZ, cioè correlaz tra blocco(colonna, lungo NxN=N^2) e se stesso (NxN^2). Ogni foglio è l'autorrelaz di
%un blocco per se stesso, quindi sono tante quante num_block
r_blocks = zeros(N^2,N^2,num_block); %num_block lo avevo già calc, ed era il num dei blocchi ceh entrano nella imm
for block_idx = 1:num_block %scandisco per layer che sono di numero uguale al num di blocchi della imm
    block = imblocks_zm(:,block_idx); %prendo il blocco, che sono le colonne di imblocks_zm infatti
    r_blocks(:,:,block_idx) = block * block'; %faccio l'autocorrelaz
    %r_blocks è la matrice che contiene ,per ogni layer, la correlazione di ogni blocco con se stesso (aka l'autocorrelaz)
end

R = mean(r_blocks,3);
% mean(A,3) is a MATRIX containing the mean of each elements of each layer (che sono i blocchi); cioè medio per fogli (blocchi) sovrapposti

%% 6) Perform KLT coding - estimate the PSNR and display the reconstructed image

% Compute transform matrix from correlation
[V,D] = eig(R); %diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors
T_klt = V; %just for convenience, dovrebbe essere V', ma per far funzionare la molt di dopo mi ritrovo solo così con le dimensioni (343)

Q = 23;
Qvector = Q*ones(N*N,1); %xke ones ???????? forse perchè voglio solo cambiare le basi di rappresentazioni della mia imm

klt_symbols = zeros(K^2,num_block); % sono le REPRODUCTION REGION, distribuiti per colonne, ogni colonna è riferita ad un blocco
block_idx = 1;

y_klt = zeros(size(Y));
% For each block, ora riprendo il blocco quadrato
for r = 1:size(Y,1)/8
    for c = 1:size(Y,2)/8
        % 6a) consider block of size 8x8
        block_r_idxs = (r-1)*N+1:r*N;
        block_c_idxs = (c-1)*N+1:c*N;

        block = Y(block_r_idxs,block_c_idxs);

        % 6b) compute the KLT
        this_block_mean = mean(block(:)); %faccio la media di tutti i valori della matrice block 8x8
        block_klt = T_klt'*(block(:)- this_block_mean); %applico la matrice di trasf al blocco unbiased, prima era dct2, e toglievo 127

        % 6c) threshold quantization, cioè QUANTIZZO
        block_klt_coeff = round(block_klt./Qvector);
        block_klt_q = block_klt_coeff.*Qvector;

        klt_symbols(:,block_idx) = block_klt_coeff; %fillo la matrice definita sopra
        block_idx = block_idx + 1;

        % 6d) reconstruct the image from quantized coefficients (use reshape() if needed)
        block_iklt = T_klt * block_klt_q + this_block_mean; %block_iklt è una 64*1 quindi uso il reshape
        %sotto per ricostruire il blocchetto che mi serve per costruire l'immagine
        y_klt(block_r_idxs,block_c_idxs) = reshape(block_iklt,8,8); %

    end
end

% display image and compute PSNR  and entropy
%
% figure()
% subplot(1,2,1);
% imshow(uint8(Y));
% title('Original');
%
% subplot(1,2,2);
% imshow(uint8(y_klt));
% title('KLT');

mse_klt = mean((Y(:)-y_klt(:)).^2);
PSNR_klt = pow2db(255.^2/mse_klt);

disp(['PSNR KLT: ' num2str(PSNR_klt)])

klt_symbols = klt_symbols(:); %unfoldo la matrice dei simboli in un vett colonna totale
klt_symbols_values = unique(klt_symbols);
klt_symbols_count = hist(klt_symbols, klt_symbols_values);
klt_symbols_prob = klt_symbols_count./sum(klt_symbols_count);
klt_entropy = -sum(klt_symbols_prob.*log2(klt_symbols_prob));

fprintf('Entropy KLT: %.4f bit/symbol\n',klt_entropy);
