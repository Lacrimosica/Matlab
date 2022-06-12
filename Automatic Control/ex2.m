s= tf('s')

H= (s+5)/(s^2 + 5*s -4)

zeros_H = zero(H)
poles_H = pole(H)

zpk(H)

% another question

s = tf('s')

U=2/s^2
H= (2*s +1)/ (s+4)^2

Y= zpk(minreal((H*U) , 1e-3))

[num, den] = tfdata(Y, 'v')

[r, p] = residue(num,den)