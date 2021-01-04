% Es 10

% 10.a
% Given a signal x(n) = [3, 2, 1, 0, 1], n in [-2, 2]
% Given a LTI system with h(n) = [1, 3, 2.5, 4, 2], n in [0, 4]
% Compute the output of the system using conv. Which is the support of y(n)?
% Exploiting the convolution theorem, compute Y(z) = X(z) H(z)
% Which is the order of polynomial H(z)?
% 10.b
% Compute the roots of H(z). 
% Write y1(n) as the convolution of x(n) with the filter cascade h(n) =
% h1(n) * h2(n) * h3(n) etc..


close all
clearvars
clc

%% Es 10.a

x = [3, 2, 1, 0, 1];
n_x = -2:2;

h = [1. 3. 2.5, 4, 2];
n_h = 0:4;

%% compute convolution using conv

y = conv(x, h);
% support of y:
n_y = n_x(1) + n_h(1):n_x(end) + n_h(end);

%% plot the signal

figure;
stem(n_y, y);

%% expression of H(z)

H_z = h;

%% order of the polynomial H(z)

order = n_h(end);

%% Es 10.b : roots of the filter

h_roots = roots(h);

%% compute y1 as the convolution of x(n) with the filter cascade.

h_0 = h(n_h == 0);

% initialize the filter cascade as a delta. 
h_cascade = 1;
% support of the delta function
sup_cascade = 0;

% loop over the roots
for r = 1:length(h_roots)
    
    % create the elementary sequence: its support is n = [0, 1]
    h_r = [1, -h_roots(r)];
    % convolve by the cascade
    h_cascade = conv(h_r, h_cascade);
    % new support of the cascade
    sup_cascade = sup_cascade(1) + 0:sup_cascade(end) + 1;
       
end

% final cascade
h_cascade = h_0 * h_cascade;

y1 = conv(x, h_cascade);
sup_y = n_x(1) + sup_cascade(1):n_x(end) + sup_cascade(end);


%% plot the result

figure;
stem(n_y, y);
hold on
stem(sup_y, real(y1), '--');
grid



