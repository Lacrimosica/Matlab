clear all
close all
clc



N=input("What is the value of N?");

n=0:N-1;

UF= zeros(N,N);

for k=0:N-1
    uk = exp(1j*2*pi*k*n/N) / sqrt(N);
    UF(:,k+1)= uk; %%indexing starts from 1
end



 figure()
 imagesc(real(UF'*UF)), title("real(UF'*UF)")
 colorbar


 figure()
 imagesc(imag(UF'*UF)), title("imag(UF'*UF)")
 colorbar

  figure()
 imagesc(real(UF*UF')), title("real(UF*UF')")
 colorbar

 figure()
 imagesc(imag(UF*UF')), title("imag(UF*UF')")
 colorbar

 %%................
 x = randn([N,1]);
 c= UF' *x
    
 n2c = norm(c,2); %%to get the square norm
 n2x = norm(x,2) %%square norm of x
 %%this is the proof of parsevals inequality

 %%-------------------


 
