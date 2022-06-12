clear all
close all
clc

s = tf('s')
TOL = 1e-3;

H= zpk(minreal( 1/( (s^2 +9) * (s^2 +2 *s +10) ), TOL))

sys=ss(H)



% supposing U is an impules , aka dirac delta
U=1;


Y= H*U;


[num, den] = tfdata(Y, 'v')

[res, poles] = residue(num, den)


% theese are the components for the inverse laplace transfrom 
% of the output response equation(i.e in time)

R = res(1);
sigma = real(poles(1))
omega = imag(poles(1))
twotheta = 2*abs(R)
psi = angle(R)



R = res(3);
sigma = real(poles(3))
omega = imag(poles(3))
twotheta = 2*abs(R)
psi = angle(R)
