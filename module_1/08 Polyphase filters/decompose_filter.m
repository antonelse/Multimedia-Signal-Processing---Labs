function [ e ] = decompose_filter(h, M )
%DECOMPOSE Polyphase decomposition of the filter h in M subfilters
    L=length(h);
    K=ceil(L/M);
    
    if K*M > L
        h(K*M)=0; % zero-pad up to K*M     
    end
    e=reshape(h,M,length(h)/M).';
    
end

