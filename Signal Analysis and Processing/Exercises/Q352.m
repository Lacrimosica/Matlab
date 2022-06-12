clear all; close all;   % initial clean up
dt=0.1;                 % time-step for numerical integration (s)
t=0:dt:160;             % time-array for numerical integration (s)

a=10+t/8;                           % acceleration (m/s^2) vs time (s)
figure;plot(t,a)    ;
                                % plotting acceleration (m/s^2) vs time (s)
title('acceleration vs time');
ylabel('m/s^2');xlabel('t');grid on;
pause;
 % waiting for keystroke to go on
%
for n=1:numel(t)                                % numerical integration to obtain speed
    if n==1
        v(1)=0;
    else               % trapezoidal integration rule
        v(n)=v(n-1)+1/2*(a(n-1)+a(n))*dt;   
    end
end


figure;plot(t,v);hold on;
 % plotting speed (m/s) vs time (s)
title('speed vs time');
ylabel('m/s');xlabel('t');grid on;
pause;
 % waiting for keystroke to go on
%
va=10*t+t.^2/16;
 % analytical formula for speed
plot(t,va,'r--','LineWidth',2);hold off; % plotting it too, for comparison,
% red dashed
pause;
 % waiting for keystroke to go on
%
for n=1:numel(t)
 % numerical integration to obtain height
if n==1
h(1)=0;
else
 % trapezoidal integration rule
h(n)=h(n-1)+1/2*(v(n-1)+v(n))*dt;
end;
end;
figure;plot(t,h);hold on;
 % plotting height (m) vs time (s)
title('altitude vs time');
ylabel('m');xlabel('t');grid on;
pause;
 % waiting for keystroke to go on
%
ha=5*t.^2+t.^3/48;
 % analytical formula for height
plot(t,ha,'r--','LineWidth',2);hold off; % plotting it too, for comparison,