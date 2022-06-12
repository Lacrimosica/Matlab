clear all
close all


%% This section devoted to a generic raised-cosine pulse

T=pi-1;
M=1001;
t=linspace(-1/2,1/2,M).*T;
roll_off=0.3;

tau=1;
y_ref=root_raised_cosine(t,roll_off,tau,0,0);

y=tau/T*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
hold on;
xlabel('time')
axis([min(t),max(t),-0.2,1.2])
grid on;
hold on;
pause;


plot(t,y,'b','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
xlabel('time')
title(['N=0']);
hold off;
figure(h);
pause;
y_ref_energy=sum(abs(y_ref).^2)*T/M;

N=100;
for n=1:N
    
    yn=tau/sqrt(T)*sinc(n*tau/T)*(cos(pi*roll_off*n*tau/T)./(1-(2*roll_off*n*tau/T)^2));
    y=yn*(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))/sqrt(T)+y;
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    axis([min(t),max(t),-0.2,1.2])
    grid on;
    %title(['N='  num2str(n) ]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);

    figure(h);
    hold off;
    
    pause;
    
end;    
pause;

%% This section dealing with the triangular signal

figure;

tau=0.5;

for i=1:1000
    y_ref(i)=triangle(0,tau,t(i));
end
y=tau/T*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
xlabel('time')
hold on;
pause;

plot(t,y,'b','LineWidth',2);
grid on;
figure(h);
hold off;
pause;
y_ref_energy=sum(abs(y_ref).^2)*T/M;

N=100;
for n=1:N
    
    y=(sin(pi*n*tau/T)/pi/n)^2*T/tau*(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    %title(['N='  num2str(n)]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
    figure(h);
    hold off;
    
    pause;
    
end;   
pause;

%% This section dealing with the semi-circular signal

figure;
tau=0.5;

for i=1:M
    if abs(t(i))<=tau
        y_ref(i)=sqrt(1-(t(i)./tau).^2);
    else
        y_ref(i)=0;
    end
end


y=pi/2*tau/T*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
xlabel('time')
hold on;
pause;

plot(t,y,'b','LineWidth',2);
grid on;
figure(h);
pause;
hold off;
y_ref_energy=sum(abs(y_ref).^2)*T/M;

N=100;
for n=1:N
    
    y=besselj(1,2*pi*n/T*tau)/2/n*(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    %title(['N='  num2str(n)]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
    figure(h);
    pause;
    hold off;    
    
end;   
pause;


%% This section of the code to represent the unilateral exponential signal

%
T=pi-1;
f0=1/T;
M=1001;

t=linspace(-1/2,1/2,M).*T;

a_const=2;

for i=1:M
    y_ref(i)=exp(-a_const*t(i)).*u(t(i));
end
y=1/sqrt(T)*1/sqrt(T)*(1-exp(-a_const*T/2))./(a_const)*ones(1,length(t));
%y=0.2*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
grid on;
xlabel('time')
hold on;
pause;


plot(t,y,'b','LineWidth',2);
grid on;
figure(h);
hold off;
pause;
y_ref_energy=sum(abs(y_ref).^2)*T/M;
N=300;
for n=1:N
    
    spn=1/sqrt(T)*(1-exp(-a_const*T/2).*((-1).^n))./(a_const+1i*2*pi*n*f0);
    smn=1/sqrt(T)*(1-exp(-a_const*T/2).*((-1).^(-n)))./(a_const-1i*2*pi*n*f0);
    y=spn/sqrt(T)*exp(j*2*pi*n/T*t)+smn/sqrt(T)*exp(-j*2*pi*n/T*t)+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    %title(['N='  num2str(n)]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
    figure(h);
    hold off;
    
    pause;
    
end;    
pause;


%% This section of the code to represent the rectangular signal


%
T=pi-1;
M=1001;
t=linspace(-1/2,1/2,M).*T;

tau=1;

for i=1:M
    y_ref(i)=HPi(0,tau,t(i));
end
y=tau/T*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
grid on;
xlabel('time')
hold on;
pause;

plot(t,y,'b','LineWidth',2);
grid on;
figure(h);
hold off;
pause;
y_ref_energy=sum(abs(y_ref).^2)*T/M;

N=300;
for n=1:N
    
    y=sin(pi*n*tau/T)/pi/n*(exp(j*2*pi*n/T*t)+exp(-j*2*pi*n/T*t))+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    %title(['N='  num2str(n)]);
    y_err=y-y_ref;
    err_energy=sum(abs(y_err).^2)*T/M;
    title(['N='  num2str(n) ' ; error energy perc '   num2str((err_energy)./y_ref_energy*100)]);
    figure(h);
    hold off;
    
    pause;
    
end;    
pause;


%% This section dealing with a delayed rectangular signal
%

figure;
tau=1;
% value of the delay
t0=0.3;

for i=1:M
    y_ref(i)=HPi(t0,tau,t(i));
%    y_ref(i)=y_ref(i)+porta(t0-T,tau,t(i));
end
y=tau/T*ones(1,length(t));
h=figure;
plot(t,y_ref,'r','LineWidth',2)
axis([min(t),max(t),-0.2,1.2])
xlabel('time')
hold on;
pause;

plot(t,y,'b','LineWidth',2);
grid on;
figure(h);
hold off;
pause;
y_ref_energy=sum(abs(y_ref).^2)*T/M;

N=100;
for n=1:N
    
    y=sin(pi*n*tau/T)/pi/n*(exp(j*2*pi*n/T*t)*exp(-j*2*pi*n/T*t0)+exp(-j*2*pi*n/T*t)*exp(j*2*pi*n/T*t0))+y;
    
    plot(t,y_ref,'r','LineWidth',2)
    axis([min(t),max(t),-0.2,1.2])
    xlabel('time')
    hold on;
    plot(t,y,'b','LineWidth',2);
    grid on;
    title(['N='  num2str(n)]);
    figure(h);
    hold off;
    
    pause;
    
end; 