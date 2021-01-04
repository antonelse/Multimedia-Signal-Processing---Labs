function [z_out, p_out, b_out, a_out] = allpass(b,a)

% to find the all-pass filter, compute first the zeros and poles of Hz
% NB: never use the name 'zeros' because it is the name of a MATLAB function!
zeroes = roots(b);
poles = roots(a);
% NB: zeros and poles are extracted in column vectors

z_out = zeroes;
p_out = poles;

% loop over poles
for p = poles
    
    % check zeros
    if any(zeroes-1/conj(p)).^2 <= 1e-4
        % if there is already one zero which is the reciprocal conjugate of the
        % pole, go on to check the next pole
        continue
    else
        % otherwise, add one zero = reciprocal conjugate of the pole
        % remember that zeros and poles are column vectors
        z_out = [z_out; 1/conj(p)];
    end
    
end

% loop over zeros
for z = zeroes
    
    % check poles
    if any(poles-1/conj(z)).^2 <= 1e-4
        % if there is already one poke which is the reciprocal conjugate of the
        % zero, go on to check the next zero
        continue
    else
        % otherwise, add one pole = reciprocal conjugate of the zero
        % remember that zeros and poles are column vectors
        p_out = [p_out; 1/conj(z)];
    end
end

% create the related polynomials in z^-1 
b_out = poly(z_out);
a_out = poly(p_out);

end

