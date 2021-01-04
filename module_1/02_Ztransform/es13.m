% Es 13

% Given y(n) = -2 y(n-1) -  y(n-2) + x(n) + 2rho cos(theta)x(n-1) +
% rho^2 x(n-2)
% rho = 0.9, theta = pi/8
%    
% Which is the expression of H(z)?
% Compute its zeros and poles.
% Plot its zeros and poles.
% Is the system stable?

close all
clearvars
clc

%% parameters

rho = 0.9;
theta = pi/8;

%% define H(z)

% denominator coefficients (from a_0 to a_D)
A_z = [1, 2, 1];
% numerator coefficients (from b_0 to b_N)
B_z = [1, 2*rho*cos(theta), rho^2];

%% find the expression of h(n)

% from theory:
n = 0:200;
h = zeros(size(n));
h(1) = rho^2;
h = h + (2*rho*cos(theta) - 2*rho^2)*(-1).^n + (1 + rho^2 - ...
    2*rho*cos(theta))*(n+1).*(-1).^n;

%%% plot the theoretical filter.
figure; stem(h);

%%% use filter.m
delta = zeros(size(n));
delta(1) = 1;
h_filter = filter(B_z, A_z, delta);

hold on;
stem(h_filter);

%%% use residuez.m
[r, p, c] = residuez(B_z, A_z);
% NB: c is ordered from c_0 to c_(N-D).
% NB: each residue is associated to one pole.
% if one pole has multiplicity > 1, the residues are ordered as m=1, m=2...

h_partial = zeros(1, length(n));
% fil the first samples of the filter with the coefficients of vector c.
h_partial(1:length(c)) = c;

for r_i = 1:length(r)
    
    if r_i > 1 && (p(r_i) - p(r_i-1)) < 1e-4        
        % elementary sequence associated to the residue
        h_el = r(r_i) *(n+1).*p(r_i).^n;
        
    else
        % elementary sequence associated to the residue
        h_el = r(r_i) * p(r_i).^n;
    end
    
    % update h_c
    h_partial = h_partial + h_el;
    
end

hold on; 
stem(h_partial, '--')

%% zeros and poles

zeroes = roots(B_z);
poles = roots(A_z);

% theoretical zeros and poles
poles_th = [-1; -1];
zeros_th = [rho * exp(1i*theta); rho*exp(-1i*theta)];

%% plot zeros and poles

figure;
zplane(zeroes, poles);



