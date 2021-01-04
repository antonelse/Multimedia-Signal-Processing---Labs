function [ x_m ] = decompose_signal( x,M )
%DECOMPOSE_SIGNAL decompose x into M components for
% polyphase filtering
    Nm=floor(length(x)/M);
    x_m=zeros(Nm,M);
    x_0=x(1:M:end);
    x_m(1:length(x_0),1)=x_0;
    for m=1:M-1
        x_=x(M-(m-1):M:end);
        x_m(2:length(x_)+1,m+1)=x_;
    end
end

