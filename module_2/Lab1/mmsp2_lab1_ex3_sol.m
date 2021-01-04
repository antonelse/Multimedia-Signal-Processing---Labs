%% MMSP2 - Lab 1
%  Exercise 3 - Discrete memoryless source coding

clearvars
close all
clc

%% 1) Generate one realization of length 1000000 of the following process:
%%    y(n)=min(max(0,round(rho*y(n-1)+w(n))),15)
%%    where rho=0.95 and w(n) is Gaussian with variance=1
%%M = max(A,[],dim) returns the maximum element along dimension dim.
%%For example, if A is a matrix, then max(A,[],2) is a column vector containing the maximum value of each row.
%MA MI DICONO CHE 15 in questo caso indichi il val massimo assumibile dal processo, quindi è una specie di clipping,
%dato che se ho una matrice in ingresso indica la dim di azione, ma se ho un vett mi indica il valore massimo
%assumibile, quindi clipping


rng(21); %credo che serva in modo tale da avere stessa sequenza a tutti gli studenti,
%una specie di randomicità controllata e prevedibile
%%rng (seed) filla il generatore di numeri casuali usando il
%%seme intero non negativo in modo che rand, randi e randn producano una sequenza prevedibile di numeri.
N = 1e6;
rho = 0.95;
w_std = 1; %%è la radice della varianza ma è sempre 1

z = randn(N,1)*w_std;

% we have two different ways to find x(n) (use only one of them):

% 1)first way: filter function, cioè riscrivo y nel dominio zeta e calcolo la fdt
A =[1, -rho];
B = 1;
yy = filter(B, A, z);
yy = min(max(0,round(yy)),15);

% second way: for loop
y = zeros(N,1);
y(1) = z(1); %%perchè potendo agire il feedback, setto la condiz iniziale
for n = 2:N
    y(n) = rho*y(n-1) + z(n);
end
y = min(max(0,round(y)), 15);

figure();
stem(y);
title('y');

%% 2) Determine the size of the alphabet of the source
%%    hint: use the function unique()

alphabet_y = unique(y);
alphabet_len = length(alphabet_y);

disp('Alphabet');
fprintf('length(alphabet): %d\n',alphabet_len);

%% 3) Find H(Y) assuming that x is a discrete memoryless source

d_y = hist(y,alphabet_y);
p_y = d_y/sum(d_y);

H_y = -sum(p_y(d_y > 0) .* log2(p_y(d_y > 0)));
fprintf('entropy of x: %.3f bit/symbol\n',H_y);

%% 4) Let J=y(n) and K=rho*y(n-1). Compute p(J,K) and H(J,K)
j = y;
k = round(rho*[0; y(1:end-1)]); %%moltiplicando per il vett colonna,
%sto costruendo il processo (IN COLONNA) escludendo y(n) (sample in 0) e andando avanti fino a n-1

alphabet_k = unique(k);

d_joint = hist3([j, k],{alphabet_y,alphabet_k});
disp(['size(count_joint): ' mat2str(size(d_joint))]);

p_joint = d_joint/sum(d_joint(:));

figure();
imagesc(db(p_joint));
title('Joint PMF');

% compute joint entropy
H_joint = -sum(sum(p_joint(d_joint > 0) .* log2(p_joint(d_joint > 0))));
fprintf('joint entropy: %.3f bit/2 symbols <= 8 bit/2 symbols\n',H_joint);


%% 5) Compute the conditional entropy H(J|K)

d_k = hist(k,alphabet_k);
p_k = d_k/sum(d_k);

H_k = -sum(p_k(d_k > 0) .* log2(p_k(d_k > 0)));
fprintf('entropy of y: %.3f bit/symbol\n',H_k);

% compute conditional entropy using chain rule
H_cond_cr = H_joint - H_k;
fprintf('cond entropy X|Y: %.3f bit/symbol\n',H_cond_cr);
