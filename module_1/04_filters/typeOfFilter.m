function [filter_type] = typeOfFilter(b,a)

% compute zeroes and poles
zeroes = roots(b);
poles = roots(a);

% stability is related to poles only
if any(abs(poles) > 1)
    filter_type = -1;
else
    % the system is stable
    % check if it is minimum phase 
    % also zeros must be inside the unit circle
    if any(abs(zeroes) > 1)
        % the system is not minimum phase
        filter_type = 0;
    else
        % the system is minimum phase
        filter_type = 1;
    end
end

end

