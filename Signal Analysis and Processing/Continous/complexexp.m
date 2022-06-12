clear all
close all
clc

t=0:0.01:10;
z=exp(j*2*pi*t);
figure();
plot3(t,real(z), imag(z),  'LineWidth', 2);
grid on;
axis equal; 
axis([0,10,-1.2,1.2, -1.2,1.2]);
xlabel('time(s)') ;
ylabel('real part');
zlabel('imaginary part'); 
hold on;
plot3(t, zeros(size(real(z))), zeros(size(imag(z))) ,  'LineWidth' , 1);
hold off;