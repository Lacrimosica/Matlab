clear all
close all
clc

s= tf('s'); 


G= 2122.4/(s*(s+59.24));
C= 4.54*(1+s/59.2)/(1+s/218.8);


t_stop = 1;

% to automatically 
% run the simulation
x = sim("mysim.slx")

figure(1)
plot(x.y.time, x.y.Data , 'b' , 'LineWidth',2)
grid on

figure(2)
plot(x.e.time , x.e.Data , '--b' ,'LineWidth', 2  )
grid on


