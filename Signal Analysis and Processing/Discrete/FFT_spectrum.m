clear all
close all
t0=1;%normalized value, t0=1 s 
Deltat=t0/4;
N=32;
t=[-N/2:N/2-1]*Deltat;
x=exp(-abs(t)/t0);
N=length(x);
X1=fft(x)*Deltat;
X=fftshift(X1);
f=[-N/2:N/2-1]/(N*Deltat);
plot(f,abs(X),'-ro',f,2*t0./(1+(2*pi*f*t0).^2),'-bo')
legend('Estimated','True')
grid on, xlabel('f (Hz)'), ylabel('|X(f)|')
