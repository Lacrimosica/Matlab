clear all
close all
clc

rho=0.99

theta_deg = 185 ; theta = theta_deg/180*pi;

% poles
p1 = rho*exp(j*theta);
p2=conj(p1)

% sampling interval
Deltat = 1 ; 

f = [-0.5:0.01:0.5]/Deltat

z= exp(j*2*pi*f*Deltat); 


H= z.^2./((z-p1).*(z-p2));

figure()
plot(f,abs(H)) ; grid on

xlabel('frequency 1/delta_t')
ylabel('|H|')

X=ones(20,1)
X(17:20)=0

Y = X.
figure()
plot(f,Y) ; grid on
