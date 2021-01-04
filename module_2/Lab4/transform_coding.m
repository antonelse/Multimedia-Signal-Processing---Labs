function [R_opt,snr] = transform_coding(x,T,R)
%TRANSFORM_CODING Folding, projection, quantization, inverse projection,
%unfolding with Optimal Bit Allocation
%   x: 1D input signal
%   T: transformation matrix (square matrix, each row is a basis)
%   R: bits per symbols budget

% Determine group of symbols length based on transformation length
N = size(T,1);
%la dim di T, e quindi N, dipende dal numero di simboli che voglio encodare a gruppi; cmq sia è una matr quadrata sempre

% reshaping (possibly cutting) the signal in order to have groups of symbols
x = x(1:floor(length(x)/N)*N); %prendo i valori che VERAMENTE mi entrano se raggruppo a gruppi di N ! scarto l'eccesso
%MA IL *N ha senso poichè mi dice il num tot di elementi alla fine della moltiplicaz. Taglio fuori gli elementi che nn riescono a formare un gruppo
X = reshape(x, N, length(x)/N); %dispone il vettore x (aka ingresso) per "num di simboli per gruppo" righe e "numero gruppi" colonne,
%così ad ogni colonna corrisponde un gruppo di simboli, dato che devo fare moltiplicazione riga di T (base) per colonna di X (gruppo di simboli)



% apply the transformation
Y = T*X;
%nb T è NxN; X è NxNumGruppi;  Y è NxNumGruppi
%la prima riga di Y indica quanto la prima base della trasformata contribuisce alla "sintesi" dei gruppi di simboli
%(il num delle colonne è il numero dei gruppi di simboli)


% compute optimal bit allocation
var_Y = var(Y,0,2); % calcolo la varianza della trasformazione per ogni colonna (2) (cioè la valuto per ogni gruppo)
%settando la normalizzazione al valore default (0)

R_opt = round(R+0.5*log2(var_Y./geomean(var_Y))); %geomean è la media geometrica
%formula di pag 55

% compute delta for each coefficient
%faccio il max e min per ogni sample di ogni singolo gruppo poichè:
%max(A,[],2) is a column vector containing the maximum value of each row.
%e la Y contiene (in colonne) i contributi di ogni singola base della trasformata per ogni gruppo
%il vett delta quindi contiene il delta associato ad ogni colonna quindi ad ogni gruppo trasformato
delta = (max(Y,[],2)-min(Y,[],2))./(2.^R_opt-1);

% repeat delta to fit the size of Y
% per fare l'operazione matriciale dopo, ricopio (replico) questo vett colonna delta in più colonne
%(uguali al num di gruppi) ottenendo una matrice
%il secondo indice di repmat mi dice che la ripetizione di x (vett colonna) è solo una rispetto alle righe,
% quindi non mette x sotto x (esempio caso (..,2,..)
%ma lo ripete per colonne, tante volte quante il num di gruppi (aka length(x)N)
%exp sta per expanded
delta_exp = repmat(delta,1,length(x)/N);

% quantize Y (è una divisione element-wise)
Y_tilde = floor(Y./delta_exp).*delta_exp + delta_exp/2;

% invert the transformation
X_tilde = T'*Y_tilde;
%matr trasf è  matrice ortogonale, quindi inversa = trasp

x_tilde = X_tilde(:); %stendo la matrice in un vettore

% compute SNR
mse = var(x-x_tilde);
snr = pow2db(var(x)/mse);

end
