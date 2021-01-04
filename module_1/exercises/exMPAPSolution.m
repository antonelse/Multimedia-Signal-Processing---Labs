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
% your code and answer here
figure;
N=512;
plot(linspace(0,1,N),...
    20*log10(abs(freqz(b,a,N))));
xlabel('Norm. \omega');
ylabel('|H(z)| [dB]');
title('Magnitude of the freq. response');
figure;
zplane(b,a);
title('Zplane')
fprintf(['H(z) is a high-pass filter, it is stable \n',...
      '(all poles in the unit circle), but not \n',...
      'invertible (zeros outside the unitary circle)\n']);
%% Listen to x
 % your code here
 sound(x,Fs);
 pause(length(x)/Fs+1);
%% Listen to y
 % your code here
 
sound(y,Fs); 
pause(length(y)/Fs+1);
%% We want to restore x from y
% What happen if we just invert b and a? Plot the result

x_hat=filter(a,b,y);
plot(x_hat);
title('x=h^{-1} y');
disp('The filter is unstable!');

%% Find a way to restore x from y: pseudo-code below

% find the minimum-phase all-pass decomposition of the filter
% z=...
% p=...
% z_mp= ...
% p_mp=...
% z_ap= ...
% p_ap=...
% b_ap=...
% a_ap=...
% b_mp=...
% a_ap=...

z=roots(b);
p=roots(a);
z_max=z(abs(z)>=1);
z_min=z(abs(z)<1);

z_ap=z_max;
p_ap=1./conj(z_max);
p_mp=p;
z_mp=[z_min; p_ap];
b_ap=poly(z_ap);
a_ap=poly(p_ap);
b_mp=poly(z_mp);
a_mp=poly(p_mp);

% compute the inverse and listen to it
% hint: H^-1=H_mp^-1 * H_ap

z_inv=[z_ap; p_mp];
b_inv=poly(z_inv);
p_inv=[p_ap; z_mp];
a_inv=poly(p_inv);

x_hat=filter(b_inv, a_inv, y);
sound(x_hat, Fs);
%% What kind of filter you are applying?

% your code and answer here
figure;
N=512;
%plot(linspace(0,1,N),...
    %20*log10(abs(freqz(b_inv,a_inv,N))));

plot(linspace(0,1,N),...
    20*log10(abs(freqz(a_mp,b_mp,N))));
hold on
plot(linspace(0,1,N),...
    20*log10(abs(freqz(b_ap,a_ap,N))));
plot(linspace(0,1,N),...
    20*log10(abs(freqz(b_inv,a_inv,N))));
legend('MP component','AP component', '|H^{-1}(z)|');

xlabel('Norm. \omega');
ylabel('|H(z)| [dB]');
title('Magnitude of the freq. response');
figure;
hold on

scatter(real(z_mp), imag(z_mp),'rx')
scatter(real(p_mp), imag(p_mp),'ro')
%zplane(a_mp,b_mp,'k');

scatter(real(p_ap), imag(p_ap),'kx')
scatter(real(z_ap), imag(z_ap),'ko')
legend('Poles of H_{mp}^{-1}',...
       'Zeros of H_{mp}^{-1}',...
       'Zeros of H_{ap}',...    
       'Poles of H_{ap}');
omega=linspace(0,2*pi,2048);
plot([0,0],[-2,2],'b--');
plot([-2,2],[0,0],'b--');
plot(cos(omega),sin(omega),'b--');
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
xlabel('Real(z)');
ylabel('Imag(z)');
hold off

title('Zplane')
fprintf(['H(z) is a low-pass filter\n']);