% Es 4

% Generate the signal x(n) = (0.8)^n u(n), n = 1:20
% Generate the signal y1(n) = x(n-5),  n = 1:20
% Generate the signal y2(n) = x(n+5), n = 1:20
% Hint: Consider using circshift instead of for loops.
% Plot the signals in the same figure.

close all
clearvars
clc

%% generate the signal

n = 1:20;
% NB: you have to put .^
x = (.8).^n;

%% shift the signal

y1 = zeros(size(x));
y2 = zeros(size(x));

% for loop
for i = 1:20
    
    if i > 5
        y1(i) = x(i-5);
    end    
    
    if i <= 15
        y2(i) = x(i+5);
    end
    
end

% you can use also 'circshift' but be careful to introduce zeros
y1_c = circshift(x, 5);
y1_c(1:5) = 0;

y2_c = circshift(x, -5);
y2_c(end-5:end) = 0;


%% plot the 3 signals in the same figure.

ll = {};
figure;
stem(x)
ll{1} = 'x(n)';
hold on, 
stem(y1)
ll{2} = 'x(n - 5)';
hold on, 
stem(y2)
ll{3} = 'x(n + 5)';
l = legend(ll, 'fontsize', 14);
grid

ll = {};
figure;
stem(x);
ll{1} = 'x(n)';
hold on,
stem(y1_c); 
ll{2} = 'x(n-5)';
hold on,
stem(y2_c); 
ll{3} = 'x(n+5)';
l = legend(ll, 'fontsize', 14);
grid



