function X=myFFT(x)
% assumption: column vectors
% assumption: x has a power of 2 number of samples
L=length(x);
if L>2
    xe=x(1:2:end);
    xo=x(2:2:end);
    Xe=myFFT(xe);
    Xo=myFFT(xo);
    k=[0:L-1]';
    X=[Xe;Xe]+exp(-1j*2*pi*k/L).*[Xo;Xo];
else
    X=[x(1)+x(2);x(1)-x(2)];
end
return
