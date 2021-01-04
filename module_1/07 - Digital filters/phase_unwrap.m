function [ output_args ] = phase_unwrap( H )
%PHASE_UNWRAP Summary of this function goes here
%   Detailed explanation goes here
    phase=unwrap(angle(H));
    for i=2:length(H)
        if abs(phase(i)-phase(i-1))>=pi            
        end
    end
    
    figure; plot(unwrap(phase)/pi);   
end

