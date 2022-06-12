clear all
close all
clc


% implmentation of one-pole systems
% H(z) = z/(z-p)
% y[n] = x[n] +py[n-1]

N =100;
n= [-N:N]';

x=zeros(2*N+1,1);
x = x.* (n==0);
figure()
stem(n,x) ; grid on
xlabel('n')
ylabel('x[n]')
% % system implementation
p=0.5;
y=x*0;
for k=2:length(x)
    y(k) = p*y(k-1) + x(k)
end

figure()
stem(n,y), grid on, hold on
xlabel('n')
ylabel('y[n]')
%theoretical impulse response

h=(p.^n).*(n>=0) ; 
stem(n,h,'x'), grid on , hold on
legend('measured' , 'theoretical')
