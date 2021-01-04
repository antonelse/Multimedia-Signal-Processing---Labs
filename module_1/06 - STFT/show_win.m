function [  ] = show_win(win, N, name)
%SHOW_WIN A tool to plot time and frequency domain of a window
% zero-padded to N
    M=length(win);
    
    
    
    figure;
    subplot(1,2,1);
    plot([0:M-1], win);    
    xlim([-1, M]);
    ylim([-0.1,1.1]);
    xlabel('n'); ylabel(['w_{', name,'}']);
    title([name ' in the time domain']);
    
    win(N)=0;   
    subplot(1,2,2);
    plot([0:N-1], 20*log10(abs(fft(win))));
    xlim([0, N/2+1]);
    xlabel('k'); ylabel(['|W_{', name,'}| [dB]']);
    title([name ' in the frequency domain']);
end

