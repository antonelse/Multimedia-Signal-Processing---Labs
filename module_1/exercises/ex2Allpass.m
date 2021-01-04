%% Allpass filters

% Given a fir filter B(z)=?2z^(?3)+2z^(?1)+3,
b=[3 2 0 -2];

%% 1) implement an all-pass filter equal to H(z)=B(z)/A(z)
    %a=?
    %z=roots(b);
    %p=1./conj(z);
    %a=poly(p);
     a=fliplr(b);  
   
%% 2) plot the poles and the zeros of the all-pass filter in the z-plane
    figure;
    %zplane(z,p);
    zplane(b,a);

%% 3) compute its frequency response from 0 to pi and plot its modulus and 
%       phase (two subplots), with normalized frequency
    figure;
    [H, w]=freqz(b,a,2048);
    subplot(2,1,1);
    plot(w/pi, abs(H));
    xlabel('\omega');
    ylabel('|H(z)|')
    subplot(2,1,2);
    plot(w/pi, angle(H));
    xlabel('\omega');
    ylabel('\angle H(z)')
 
%% 4) compute and plot its impulse response: is it accurate?
    h= filter(b,a,[1,zeros(1,2047)]);
    figure;
    plot([0:length(h)-1],h);
    xlabel('n');
    ylabel('h(n)');