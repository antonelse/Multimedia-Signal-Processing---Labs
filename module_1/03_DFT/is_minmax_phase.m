function [] = is_minmax_phase(zero)

if zero >= 1
    disp('Maximum phase');
elseif zero >0 && zero < 1
    disp('Minimum phase');
else
    disp('The zero is negative');
end
    
end

