% Es 8

% Given x(n) = [3, 11, 7, 0, -1, 4, 2] , n in [-3, 3]
% Given h(n) = [2, 3, 0, -5, 2, 1], n in [-1, 4]
% Compute y(n) as x(n) convolved with h(n), n in  [-7, 7]. 

close all
clearvars
clc

%% signals

% let us define the signals in the same n-domain. --> put zeros where
% signals are not defined.

n_x = -3:3;
n_h = -1:4;
x = [0, 0, 0, 0, 3, 11, 7, 0, -1, 4, 2, 0, 0, 0, 0];
h = [0, 0, 0, 0, 0, 0, 2, 3, 0, -5, 2, 1, 0, 0, 0];
n = -7:7;

%% stem the two signals

subplot(211)
stem(n, x)
title('x(n)');
grid
subplot(212)
stem(n, h)
title('h(n)');
grid

%% convolution

% fold h --> h(k) --> h(-k) (with respect to 0)
h_fold = fliplr(h);

y = zeros(size(n));

% loop over all possible n-values
cnt = 1;
for n_i = n
    
    h_shifted = circshift(h_fold, n_i);
    if n_i > 0
        h_shifted(1:n_i) = 0;
    else
        h_shifted(end + n_i :end) = 0;        
    end
    
    % scalar product between x and shift(h_fold)
    y(cnt) = x*h_shifted';
    
    % increase the bin counter for y
    cnt = cnt + 1;
    
end

figure;
stem(n, y)
title('y(n) without using conv');

%% use the matlab function 'conv'

% default:'full' convolution, with length = length(x) + length(h) - 1 
y1 = conv(x, h);
supp_full = n(1) + n(1):n(end) + n(end);

figure;
stem(supp_full, y1);
title('y(n) using full conv');

% which is the actual support of y1?
supp_true = n_x(1) + n_h(1):n_x(end) + n_h(end);

figure; 
stem(supp_true, y1(supp_full>=supp_true(1) & supp_full<=supp_true(end)));
title('y(n) inside the true support');





