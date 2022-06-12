% this simple program shows that the sum of many complex exponentials tends
% to generate a sinc function, which in turns tends to "become a delta"

clear all;
close all;

display(' ')
display('to stop execution, hit Ctrl+C')
display(' ')
display('to keep increasing B, simply press enter')

T=10;
M=1001;
t=linspace(-1/2,1/2,M).*T;

B=0;
y=zeros(1,length(t));
N=10000000;
h=figure;
plot(t,y,'linewidth',2);
axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
grid on;
figure(h);
xlabel('time [s]')
title(['B = 0']);
pause;
gcah=gca;
gcah.FontSize=16;

% B=1;
% y=B*sinc(B*t);
% plot(t,y,'linewidth',2);
% axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
% grid on;
% figure(h);
% xlabel('time [s]')
% title(['T = '  num2str(T) '   --   N = 0' ]);
% pause;

for n=1:N
    B=B+1
%   y=(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))+y;
%   equivalent form:
    y=B*sinc(B*t);

    plot(t,y,'linewidth',2);
    axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
    grid on;
    title(['B = ' num2str(B)]);
    figure(h);
    xlabel('time [s]')
    pause;
    
end;    

