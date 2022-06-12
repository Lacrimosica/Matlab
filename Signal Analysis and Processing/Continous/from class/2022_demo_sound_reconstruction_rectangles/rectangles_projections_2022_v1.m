clear all
close all

load test_signal2    % vasca, 7.5 secondi 

N=2^14;     % number of rectangles used
s=zeros(1,N);
t0=t_end/N;
f0=1/t0;
sqrt_t0=sqrt(t0);
signal= @(t) test_signal2(t)' ./sqrt_t0;


for n=1:N
if mod(n,1000)==0, display(n); end;
s(n)=integral( signal, (n-1)*t0 , n*t0, 'RelTol', 1e-6 );
end;


P_signal=@(t) (test_signal2(t)').^2 ;
signal_energy=integral( P_signal , 0, t_end, 'RelTol', 1e-11)
approximant_energy=sum(s.^2)
(signal_energy-approximant_energy)/signal_energy
10*log10((signal_energy-approximant_energy)/signal_energy)


save results2_14.mat
return