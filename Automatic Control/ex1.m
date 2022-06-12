clear all
close all
clc 

A= [0 1; -2  -3];
B=[1;0];
C=[0 1];

x0= [2 ; 2];





s = tf('s');
U=2/s;


TOL = 1e-3;
X = minreal(inv(s*eye(2)-A)*( B*U + x0 ), TOL)


% partial fraction 

[num_X1, den_X1] = tfdata(X(1) , 'v')

% tf data extracts numertor and denominator
% option v is to save them in vectors

[residues_1, poles_1] = residue(num_X1, den_X1)  

% computes residues (coeficients of partial fraction expansion
% and poles which are the roots of the denominator
% the residue function only takes values that are
% in the format that tfdata returns to you


% same thing for the socond fraction

[num_X2 , den_X2] = tfdata(X(2), 'v')

[residues_2, poles_2] = residue(num_X2, den_X2)


% each residue is corelated to the same pole in that order

A= [-3 2; -2 -3]; B=[1;0];

Y = zpk(minreal(C*inv(s*eye(2)-A)*( B*U + x0 ), TOL))


[num,den] = tfdata(Y,'v')
[residues, poles] = residue(num, den)


% for getting euler formula
theta=abs(residues(1))
psi = angle(residues(1))