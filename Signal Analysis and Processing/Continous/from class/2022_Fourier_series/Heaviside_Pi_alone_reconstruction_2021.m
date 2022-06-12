% this simple program reconstructs a Heaviside Pi by keeping on adding
% Fourier components

clear all;
close all;

display(' ')
display('to stop execution, hit Ctrl+C')
display(' ')
display('to keep on adding components, simply press enter')
display(' ')
td=input('enter the signal delay: [default 0]   ');
if isempty(td)
    td =0;
end


%%
%T=pi-1;
T0=pi;
M=1001;
t=linspace(-1/2,1/2,M).*T0;

T=1;

for i=1:M
    y_ref(i)=HPi(td,T,t(i));
end
y=T/T0*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
grid on;
xlabel('time')
hold on;
pause;

plot(t,y,'b','LineWidth',2);
grid on;
y_err=y-y_ref;
y_ref_energy=sum(abs(y_ref).^2)*T0/M;
err_energy=sum(abs(y_err).^2)*T0/M;
title(['N=1', ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
figure(h);
hold off;
pause;


N=300;
for n=1:N
    
    y=sin(pi*n*T/T0)/pi/n*(exp(j*2*pi*n/T0*t)*exp(-j*2*pi*n/T0*td)+exp(-j*2*pi*n/T0*t)*exp(j*2*pi*n/T0*td))+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    %title(['N='  num2str(n)]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T0/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
    figure(h);
    hold off;
    
    pause;
    
end;    