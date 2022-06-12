clear all
close all
clc
B=[0; 0 ];
% D=[0];
% C=[];
% sys = ss(A,B,C,D);

s=tf('s'); 
A=[ 0 1; -2 -3]
X0=[2; 2];

X=minreal(inv(s * eye(2) - A) *X0, 1e-3)

[num1 , den1] = tfdata(X(1), 'v');

[r1, p1] = residue(num1, den1)

[num2 , den2] =  tfdata(X(2), 'v')

[r2, p2] = residue(num2, den2)


