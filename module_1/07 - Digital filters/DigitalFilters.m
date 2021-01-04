clear all
clc
close all

%% Filter characterization: how to pass from one form to another

b=[3,2,1,2,3];
a=[1, -0.7, 0.2];

% Using Freqz
[H, f]=freqz(b,a);
figure(1);
subplot(2,1,1);
plot(f/pi,abs(H));
xlabel('\omega [\pi units]'); ylabel('|H|');
subplot(2,1,2);
plot(f/pi,angle(H));
xlabel('\omega [\pi units]'); ylabel('\angle H');

figure(2);
subplot(1,2,1);
zplane(b,a);
title('Difference Equation');

% Zero-pole factorization
z=roots(b);
p=roots(a);

figure(2)
subplot(1,2,2);
zplane(z,p);
title('Zero Pole Factorization');


% Unit circle estimation

H_zp=b(1)*prod((1-z*exp(-1i*f')),1)./prod((1-p*exp(-1i*f')),1);
H_zp=H_zp.';
figure(1);
subplot(2,1,1);
hold on; plot(f/pi,1+abs(H_zp.'),'r--'); %1+abs(H_zp) for better visualization
legend('Difference equation', 'Zero-Pole Factorization');
subplot(2,1,2);
hold on; plot(f/pi,0.1+angle(H_zp.'),'r--');
legend('Difference equation', 'Zero-Pole Factorization');

% Using filter and FFT
h=filter(b,a,[1 zeros(1,1000)]);
H_fft=fft(h,length(H)*2);
figure(1); 

subplot(2,1,1);
title('Transfer Function');
hold on;
plot(f/pi,abs(-1+H_fft(1:length(H))),'g--');
hold off 
legend('Difference equation', 'Zero-Pole Factorization','FFT');
subplot(2,1,2);
hold on; plot(f/pi,-0.1+angle(H_fft(1:length(H))),'g--');
hold off
legend('Difference equation', 'Zero-Pole Factorization','FFT');


% Using IFFT from freqz and estimation
figure(3);
plot(h(1:50),'b');
H=[H; flipud(conj(H(2:end)))];
h_ifft=ifft(H);
H_zp=[H_zp; flipud(conj(H_zp(2:end)))];
h_zp=ifft(H_zp);
hold on;
plot(h_ifft(1:50)+0.1,'r--');
plot(h_zp(1:50)-0.1,'g--');
hold off
xlabel('n'); ylabel('h(n)');
title('Impulse Response');
legend('Filter','IFFT from freqz','IFFT from zero-pole');
%% Filters
clear all
% one zero
show_filter([0.8],[0]);
% one pole
show_filter([0],[0.8]);

%% Low Pass
close all
%
show_filter([-1],[0.9]);

%% High Pass
close all
%
show_filter([1],[-0.9]);

%% Band Pass
close all
clear all

rho=0.9; % enhancement
phi=pi/4; % frequency


show_filter([0 0],[rho*exp(1i*phi),rho*exp(-1i*phi)]);

%%
close all
clear all

rho=0.9; % enhancement
phi=[pi/4, pi/2, pi*0.75]; % frequency


show_filter([1 -1],[rho*exp(1i*phi),rho*exp(-1i*phi)]);


%% Tuning rho
close all
for rho=0.1:0.1:1.1
    show_filter([1 -1],[rho*exp(1i*phi),rho*exp(-1i*phi)],'keep');
    pause(2);
end

%% Stopband
close all
clear all

rho=0.9; % enhancement
phi=pi/4; % frequency


show_filter([exp(1i*phi),exp(-1i*phi)],[rho*exp(1i*phi),rho*exp(-1i*phi)]);


%% Minimum phase filter
close all

%B=[1 c]
%A=[1, -a, -b];

z1=0.9; c=-z1;
p1=0.3.*exp(1i*pi/4);
p2=conj(p1);
a=p1+p2;
b=-p1*p2;

B=[1 c];
A=[1, -a, -b];
figure;
zplane(B,A);


%% Allpass filter
close all

rho=0.8; phi=0.3*pi; p1=rho*exp(1i*phi);
p=[p1, conj(p1)]; z=1./conj(p); 

show_filter(z,p);

a=poly(p);
b=fliplr(conj(a));
[H,f]=freqz(b,a);
figure; plot(f/pi,abs(H));
xlabel( 'frequency in \pi units');
ylabel( '|H(w)|')
title('Transfer Function of Allpass Filter');


%% Minphase-Allpass decomposition
clear
clc
close all

z=[2; 3*exp(1i*pi/8); 3*exp(-1i*pi/8); 0.5*exp(1i*pi/4); 0.5*exp(-1i*pi/4)]; 
p = [0.9; 0.8*exp(1i*pi/2); 0.8*exp(-1i*pi/2)];
b=poly(z);
a=poly(p);
[H,f]=freqz(b,a);

z_min=z(abs(z)<1);
z_max=z(abs(z)>=1);
p_max=1./conj(z_max);

z_mp=[z_min; p_max];
p_mp=p;
b_mp=poly(z_mp);
a_mp=poly(p_mp);
H_mp=freqz(b_mp,a);

b_ap=poly(z_max);
a_ap=poly(p_max);

H_ap=freqz(b_ap,a_ap);

H_=H_mp.*H_ap;
figure;
subplot(2,1,1);
plot(f/pi,abs(H_));
hold on;
plot(f/pi,abs(H),'--');
hold off
title('MinPhase-Allpass Decomposition');
ylabel('|H|');
xlabel( 'frequency in \pi units');
subplot(2,1,2);
plot(f/pi,angle(H_));
hold on;
plot(f/pi,angle(H),'--');
hold off
title('MinPhase-Allpass Decomposition');
ylabel('\angle H');
xlabel( 'frequency in \pi units');

