function [xq,RSB] = unifquant (x,l,Amin=min(x),Amax=max(x))
q =(Amax - Amin)/(2^l);
xq=zeros(size(x));
for n=1:size(x)(1)
  if (floor(x(n,1)/q)*q <= x(n,1)) && (x(n,1)< ((floor(x(n,1)/q) + 1)*q))
    xq(n,1) = floor(x(n,1)/q)*q + q/2;
  endif
endfor
RSB = 10*log10(mean(x.^2)/((q^2)/12));
endfunction
