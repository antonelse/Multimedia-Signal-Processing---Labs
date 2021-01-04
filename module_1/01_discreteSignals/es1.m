% Es 1

% double comment (%%) are used to divide your script
% into blocks. You can execute a block at a time by using
% "Run Section" (ctrl/cmd + enter)

% Build a signal x(n) as the sum of three different sinusoids sin(2pift) 
% at the normalized frequencies omega_1=0.11, omega_2=0.09,omega_3=0.3. 
% The sampling period is T = 0.3 seconds, and the signal is defined for 
% t in [0, 100] seconds. 

% Begin always with these three lines:
close all % close figures
clearvars % clear workspace
clc % clear command window

%% create the signal x(n)

T = .3;
duration = 100;
time = 0:T:duration;

Omega = [0.11, 0.09, 0.3];
S = zeros(length(Omega), length(time));

for o = 1:length(Omega)
    
%     for cnt = 1:length(time)
%         
%         S(o, cnt) = sin(Omega(o)*time(cnt) /T);
%         
%     end
    
    S(o,:) = sin(Omega(o) * time / T);
    
end

% % Loop can be also totally avoided (transpose Omega and multiply element-wise by time)
% S = sin(Omega'.*time./T);

x = sum(S, 1);

%% plot

figure;
plot(time, x);
xlabel('Seconds');
title('x(t=n*T)');
grid
set(gca, 'fontsize', 16)

%% period of the sinusoids

% NB: here you have to put ./ otherwise MATLAB reports an error.
P = 2*pi*T ./ Omega;

%% build x1(n)

x1 = sum(sin(2*pi*time./P'), 1);

hold on; 
plot(time, x1, '--');







