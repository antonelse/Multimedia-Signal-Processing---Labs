function [ peaks ] = findPeaks( X_abs, n )
%FINDPEAKS Find the n highest peaks in the frame X_abs
%   Very rough implementation!
    peaks=[];
    for k=1:n
        [maxX, maxI]=max(X_abs);
        peaks=[peaks, maxI];
        X_abs(maxI-5:maxI+5)=0;
    end
end

