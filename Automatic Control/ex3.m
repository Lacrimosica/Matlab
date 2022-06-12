A=[-3 2;-2 -3]; B=[1;0]; C=[0 1]; D=0;
sys = ss(A,B,C,D);



% getting the transfer function from state space representation
H=tf(sys)


% to get the state space representation from transfer function
sys2 = ss(H)

sys2.A