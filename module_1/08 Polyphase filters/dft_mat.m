function [X,w]=dft_mat(x,N)
%    if length(x)<N
%        x(N)=0;
%    end
  k=0:N-1; n=k; x=x.';
  W=exp(-1i*(2*pi/N)*(k.')*n);
  X=W*x; X=X.';
  w=k*2*pi/N;
end
