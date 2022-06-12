clear all
close all
clc

z1=-1 %first zero
zeros = [z1]; %array with all the zeros
% phi3=5
% z2

%definition of the subset og the z-plane for which H(z) is 
X= [-1.5:0.05:1.5];
Y=X;

[x,Y] = meshgrid(X,Y);
z= X+1j*y; Hz=ones(size(z));

%evaluation of H(z)
for k= 1: length(zeros)
    Hz= Hz.*(z-zeros(k));
end
for k=1:length(poles)
    Hz = Hz./(z-poles(k));
end

Hzabs = abs(Hz);
