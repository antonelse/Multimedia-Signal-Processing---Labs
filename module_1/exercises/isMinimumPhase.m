function [ result ] = isMinimumPhase( bz, ap )
%ISMINIMUMPHASE checks whether the input filter is
% a minimum phase filter
%   Detailed explanation goes here

    if length(size(bz))==2 && length(size(ap))==2
        if size(bz,1)==1 && size(bz,2)>1 &&...
           size(ap,1)==1 && size(ap,2)>1
% 1. If bz and ap are row vectors, it considers them as the coefficients 
%    of the difference equations b and a       
            b=bz; a=ap;
            z=roots(b); p=roots(a);
        elseif size(bz,2)==1 && size(bz,1)>1 &&...
           size(ap,2)==1 && size(ap,1)>1
% 2. If bz and ap are column vectors, it
%    considers them as the zeros and the poles, respectively
            z=bz; p=ap;
            b=poly(z); a=poly(p);
        else

            % 3. It returns -1 if the inputs are not correct 
            %    (input must be either two row vectors, or two column vectors)            
            result= -1;
            return
        end
    else % not row or column vectors
        result= -1;
        return
    end
% 4. It returns 1 if the filter is minimum-phase
%     and zero if it is not
    mag_z=abs(z);
    mag_p=abs(p);
    mags=[mag_z(:); mag_p(:)];
    if isempty(find(mags>1))
        result=1;
    else
        result=0;
    end


end

