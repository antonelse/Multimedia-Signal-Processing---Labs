clear all
close all
clc

%Given a filter H(z)=B(z)/A(z) write a function 
%   isMinimumPhase (bz,ap) such that:
% 1. If bz and ap are row vectors, it considers them as the coefficients 
%    of the difference equations b and a

z_mp=[0.7, 0.8, 0.6*exp(1i*pi/3),0.6*exp(-1i*pi/3)].';
z_nmp=[0.7, 0.8, 1.1*exp(1i*pi/3),1.1*exp(-1i*pi/3)].';

p_mp=[0.2, 0.5, 0.3*exp(1i*3/4*pi),0.3*exp(-1i*3/4*pi)].';
p_nmp=[0.7, 1.2, 0.3*exp(1i*3/4*pi),0.3*exp(-1i*3/4*pi)].';


% 2. If bz and ap are column vectors, it
%    considers them as the zeros and the poles, respectively

b_mp=poly(z_mp);
b_nmp=poly(z_nmp);
a_mp=poly(p_mp);
a_nmp=poly(p_nmp);

% 3. It returns -1 if the inputs are not correct 
%    (input must be either two row vectors, or two column vectors)

b_mat=randn(3,2);

fprintf('isMinimumPhase with a matrix and a vector: result is %d\n',isMinimumPhase(b_mat,a_mp));
fprintf('isMinimumPhase with a row and a column vector: result is %d\n',isMinimumPhase(b_mp,p_mp));
% 4. It returns 1 if the filter is minimum-phase
%     and zero if it is not

fprintf('isMinimumPhase with two row vectors: coefficients of a non-minimumphase filter: result is %d\n',isMinimumPhase(b_mp,a_nmp));
fprintf('isMinimumPhase with two row vectors: coefficients of a minimumphase filter: result is %d\n',isMinimumPhase(b_mp,a_mp));

fprintf('isMinimumPhase with two column vectors: coefficients of a non-minimumphase filter: result is %d\n',isMinimumPhase(z_nmp,p_mp));
fprintf('isMinimumPhase with two column vectors: coefficients of a minimumphase filter: result is %d\n',isMinimumPhase(z_mp,p_mp));