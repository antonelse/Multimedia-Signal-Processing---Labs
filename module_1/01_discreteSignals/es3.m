% Es 3 from 10/09/2018

% Generate 5 cosine tones...

close all
clearvars
clc

%% define parameters

A = [1, 0.75, 0.5, 0.25, 0.125];
Fo = 220 *[1:5];
Phi_deg = 0:45:180;
% convert the phase in radians (pi:180=Phi_rad:Phi_deg)
Phi_rad = Phi_deg * pi / 180;

duration = 1;
Fs = 8000;
time = 0:1/Fs:duration;

%% generate the five signals
% you can use for loops, or write manually every signal, or write this:

x = A'.*cos(2*pi*Fo'.*time + Phi_rad');

%% generate x6 as the sum over time of the five signals

x6 = sum(x, 1);
