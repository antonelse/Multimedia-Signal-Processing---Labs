% Es 11.b

% Given a LTI system with this transfer function: (see slides)
% Find its partial fract expansion: 
% Save in a vector r the residues
% Save in a vector p the poles
% Save in a vector c the coefficients of the remaining polynomial
% Find h(n) using filter
% Find h(n) as the cascade of filters founded with the partial fract expansion.

close all
clearvars
clc

%% define H(z)

% denominator coefficients (from a_0 to a_D)
A_z = [-1, -1, 1, 1];
% numerator coefficients (from b_0 to b_N)
B_z = [-3, -2, 1, 6, 4, 1];

%% find residues

[r, p, c] = residuez(B_z, A_z);
% NB: c is ordered from c_0 to c_(N-D).
% NB: each residue is associated to one pole.

%% find h(n) using filter.

% h(n) is the IMPULSE response to the system --> x(n) must be a delta.
n = 0:100;
delta = zeros(1, length(n));
delta(1) = 1;
h = filter(B_z, A_z, delta);

%% find h(n) as the sum of elementary filters found with the partial fract exp.

h_partial = zeros(1, length(n));
% fil the first samples of the filter with the coefficients of vector c.
h_partial(1:length(c)) = c;

for r_i = 1:length(r)
    
    if r_i > 1 && (p(r_i) - p(r_i-1))^2 < 1e-6
        % elementary sequence associated to the residue
        h_el = r(r_i) *(n+1).*p(r_i).^n;
        
    else
        % elementary sequence associated to the residue
        h_el = r(r_i) * p(r_i).^n;
    end
    
    % update h_c
    h_partial = h_partial + h_el;
    
end

figure,
stem(n, h);
hold on; 
stem(n, h_partial, '--')









