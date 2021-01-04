% Es 22

% Given the filters H1 and H2 (see slides)
% Derive A(z) and B(z) for H1 and H2 with MATLAB 
% Plot the zeros and the poles in the Z-plane using zplane
% Plot in the same figure the magnitude responses 
% as a function of normalized omega, using N = 1024 samples
% How are the magnitudes related? Why?

close all
clearvars
clc

%% define filters

B1 = conv([2, -2], [1, .5]);
A1 = conv([1, -0.8*exp(1i*pi/4)], [1, -0.8*exp(-1i*pi/4)]);

% otherwise, first derive zeros and poles and then use function
% poly which builds the related polynomial in z^-1
% NB: always consider zeros and poles in column vectors!
% z_1 = [1; -0.5];
% p_1 = [.8*exp(1i*pi/4); .8*exp(-1i*pi/4)];
% B1 = 2*poly(z_1);
% A1 = poly(p_1);

B2 = conv([1, -1], [1, 2]);
A2 = conv([1, -0.8*exp(1i*pi/4)], [1, -0.8*exp(-1i*pi/4)]);

%% plot zeros and poles in Zplane

figure;
zplane(B1, A1);
figure;
zplane(B2, A2);

%% plot the magnitude response as a function of omega

N = 1024;
[H1, omega] = freqz(B1, A1, N, 'whole');
[H2, ~] = freqz(B2, A2, N, 'whole');

figure;
plot(omega, abs(H1));
hold on; 
plot(omega, abs(H2), '--');


