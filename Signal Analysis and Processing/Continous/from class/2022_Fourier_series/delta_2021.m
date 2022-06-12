% this simple program shows that the sum of many complex exponentials tends
% to generate a sinc function, which in turns tends to "become a delta"

display(' ')
display('to stop execution, hit Ctrl+C')
display(' ')
display('to keep on adding contributions, simply press enter')

T=3;
M=1001;
t=linspace(-1/2,1/2,M).*T;

y=zeros(1,length(t));
N=10000000;
h=figure;
plot(t,y,'linewidth',2);
axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
grid on;
figure(h);
xlabel('time [s]')
title(['T = '  num2str(T) ]);
pause;

y=1/T*ones(1,length(t));
plot(t,y,'linewidth',2);
axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
grid on;
figure(h);
xlabel('time [s]')
title(['T = '  num2str(T) '   --   N = 0' ]);
pause;
for n=1:N
    
%   y=(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))+y;
%   equivalent form:
    y=2/T*cos(2*pi*n/T*t)+y;

    plot(t,y,'linewidth',2);
    axis([-T/2, T/2,min(-20,min(y)), max(50,max(y))])
    grid on;
    title(['T = '  num2str(T) '   --   N = ' num2str(n)]);
    figure(h);
    xlabel('time [s]')
    pause;
    
end;    

