clc
clear all
close all


%% Filtering x 
% We have the audio signal x, and 
% we filter through H(z)=B(z)/A(z), obtaining y
[x, Fs]=audioread('Toms_diner_16.wav');
b=[1, -1.5173, -0.0121, 0.7863, 0.1440];
a=[1. 0, 0.5 0 0.24, 0, 0.12];
y=filter(b,a,x);

%% what kind of filter is H?
% make some plots and give your answer
% your code and answer here
%% Listen to x
 % your code here
 
%% Listen to y
 % your code here
 

%% We want to restore x from y
% What happen if we just switch
% b and a? Plot the result


%% Find a way to restore x from y: pseudo-code below

% find the minimum-phase all-pass decomposition of the filter
% z=...
% p=...
% z_mp= ... --> b_mp=...
% p_mp=... -->  a_mp=...
% z_ap= ... --> b_ap=...
% p_ap=... --> a_ap=...

% compute the inverse filter
% x_hat=h^-1 * y and listen to it
% hint: H^-1=H_mp^-1 * H_ap


%b_inv=...
%a_inv=...

%x_hat=filter(b_inv, a_inv, y);

%% What kind of filter you are applying?

% your code and answer here
