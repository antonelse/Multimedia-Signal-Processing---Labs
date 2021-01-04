function [X,w]=dft(x, n)
  % modified version of the dft to use n<0
      
  N=length(n);
  k=0:N-1; 
  if size(x,1)<size(x,2)
    x=x.';   
  end
  
  W=exp(-1i*(2*pi/N)*(k.')*n);
  X=W*x; X=X.';
  w=k*2*pi/N;
end
