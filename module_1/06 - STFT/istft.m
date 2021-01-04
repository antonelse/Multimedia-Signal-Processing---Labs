function [x]= istft(X, Noverlap)
    [F, T]=size(X);
    x=zeros((T+2)*(Noverlap),1); %% just to be sure
    X_=[X;flipud(conj(X(2:end-1,:)))];
    N=(F-1)*2;
    k=1;
    for t =1:T
        x_=real(ifft(X_(:,t)));  
        x(k:k-1+N)=x(k:k-1+N)+x_;
        k=k-1+Noverlap;       
    end
    


end