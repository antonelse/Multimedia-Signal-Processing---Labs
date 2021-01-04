function [  ] = show_filter( z,p, varargin )
%SHOW_FILTER shows z-plane and transfer function of the filter
a=poly(p);
b=poly(z);
[H, f]=freqz(b,a);
f=[-flipud((f(2:end)));f];
H=[flipud(conj(H(2:end)));H];
if nargin>2 && strcmp(varargin(1),'keep')
    ;
else
    figure;
end
subplot(1,2,1);
zplane(b,a);
axis(1.2*[-1,1,-1,1])
title('zplane');
subplot(1,2,2);

[ax,p1,p2] = plotyy(f/pi,abs(H),f/pi,angle(H));
ylabel(ax(1),'|H|'); % label left y-axis
ylabel(ax(2),'\angle H'); % label right y-axis
xlabel(ax(1),'\omega [\pi units]'); % label x-axis
xlabel(ax(2),'\omega [\pi units]'); % label x-axis
title('Transfer Function');
end

