% convolution product v. 2.0
% Nov 2019

clear all;
hold off;
close all;

disp('   ');
disp('rectangular signal: ''R'' ');
disp('triangular signal: ''T'' ');
disp('unilateral decreasing exponential: ''eu'' ');
disp('power of t times unilateral decreasing exponential: ''te'' ');
disp('damped cosine: ''dc'' ');
disp('damped sine: ''ds'' ');
disp('gaussian: ''g'' ');
disp('   ');




tf1=input('which is the input signal?  ', 's');

if tf1=='R'
    t0_1= input('center time instant ? ');
    T_FWHM= input('full (time) width half maximum ?  ');
    T_1=T_FWHM; % new conventional parameter of rectangle
    f0_1= input('modulation frequency ?  ');
    c_1= input('multiplying constant ?  ');
end
if  tf1=='T' 
    t0_1= input('center time instant ? ');
    T_FWHM= input('full (time) width half maximum ?  ');
    T_1=T_FWHM; % new conventional parameter of triangle
    f0_1= input('modulation frequency ?  ');
    c_1= input('multiplying constant ?  ');
end
if tf1=='eu' 
    t0_1= input('starting time ?  ');
    a_1= input('exponent constant ''a'' ?  ');
    c_1= input('multiplying constant ?  ');
end
if tf1=='ds'
    t0_1= input('starting time ?  ');
    a_1= input('exponent constant ''a'' ?  ');
    f0_1= input('sine frequency ?  ');
    c_1= input('multiplying constant ?  ');
end
if tf1=='dc'
    t0_1= input('starting time ?  ');
    a_1= input('exponent constant ''a'' ?  ');
    f0_1= input('cosine frequency ?  ');
    c_1= input('multiplying constant ?  ');
end
if tf1=='g'
    t0_1= input('center time instant ?  ');
    a_1= input('variance ?  ');   
    f0_1= input('modulation frequency ?  ');
    c_1= input('multiplying constant ?  ');
end
if tf1=='te'
    t0_1= input('starting time ?  ');
    a_1= input('exponent constant ''a'' ?  ');
    pow_1= input('power of t ?  ');
    f0_1= input('modulation frequency ?  ');
    c_1= input('multiplying constant ?  ');
end


tf2=input('which is the impulse response signal?  ', 's');

if tf2=='R'
    t0_2= input('center time instant ? ');
    T_FWHM= input('full (time) width half maximum ?  ');
    T_2=T_FWHM; % new conventional parameter of rectangle
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end
if  tf2=='T' 
    t0_2= input('center time instant ? ');
    T_FWHM= input('full (time) width half maximum ?  ');
    T_2=T_FWHM; % new conventional parameter of triangle
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end
if tf2=='eu' 
    t0_2= input('starting time ?  ');
    a_2= input('exponent constant ''a'' ?  ');
    c_2= input('multiplying constant ?  ');
end
if tf2=='ds'
    t0_2= input('starting time ?  ');
    a_2= input('exponent constant ''a'' ?  ');
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end
if tf2=='dc'
    t0_2= input('starting time ?  ');
    a_2= input('exponent constant ''a'' ?  ');
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end
if tf2=='g'
    t0_2= input('center time instant ?  ');
    a_2= input('variance ?  ');
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end
if tf2=='te'
    t0_2= input('starting time ?  ');
    a_2= input('exponent constant ''a'' ?  ');
    pow_2= input('power of t ?  ');
    f0_2= input('modulation frequency ?  ');
    c_2= input('multiplying constant ?  ');
end




t_inf= input('lower bound of the plotted time range?  ');
t_sup= input('upper bound of the plotted time range?  ');
%N = input('numero di punti nella finestra ?  ');
N=10001;


t=linspace(t_inf, t_sup, N);
delta_t=t(2)-t(1);

switch tf1
case 'R'
    for i=1:N
        f1(i)=porta_f(t0_1,T_1,f0_1,t(i));
    end;
case 'T'
    for i=1:N
        f1(i)=triangolo_f(t0_1,T_1,f0_1,t(i));
    end
case 'eu'
    for i=1:N
       f1(i)=exp_smorzato(t0_1,a_1,t(i));
    end
case 'dc'
    for i=1:N
        f1(i)=coseno_smorzato(t0_1,a_1,f0_1,t(i));
    end;
case 'ds'
    for i=1:N
        f1(i)=seno_smorzato(t0_1,a_1,f0_1,t(i));
    end;
case 'g'
    for i=1:N
        f1(i)=gaussiana(t0_1,a_1,f0_1,t(i));
    end;
case 'te'
    for i=1:N
        f1(i)=exp_smorzato_powt(t0_1,a_1,pow_1,f0_1,t(i));
    end;
end;
% applying multiplying constant
f1=c_1*f1;
      
h1=figure(1);
set(h1,'Position',[100 100 1200, 700])
plot(t,f1,'r','linewidth',2);
title('input signal', 'fontsize',14)
xlabel('time', 'fontsize',14)
grid on;
zoom on;
figure(h1);
pause;



switch tf2
case 'R'
    for i=1:N
        f2(i)=porta_f(t0_2,T_2,f0_2,t(i));
    end;
case 'T'
    for i=1:N
        f2(i)=triangolo_f(t0_2,T_2,f0_2,t(i));
    end
case 'eu'
    for i=1:N
       f2(i)=exp_smorzato(t0_2,a_2,t(i));
    end
case 'dc'
    for i=1:N
        f2(i)=coseno_smorzato(t0_2,a_2,f0_2,t(i));
    end;
case 'ds'
    for i=1:N
        f2(i)=seno_smorzato(t0_2,a_2,f0_2,t(i));
    end;
case 'g'
    for i=1:N
        f2(i)=gaussiana(t0_2,a_2,f0_2,t(i));
    end;
case 'te'
    for i=1:N
        f2(i)=exp_smorzato_powt(t0_2,a_2,pow_2,f0_2,t(i));
    end;
end;
% applying multiplying constant
f2=c_2*f2;

h2=figure(2);
set(h2,'Position',[100 100 1200, 700])
plot(t,f2,'g','linewidth',2);
title('impulse response signal', 'fontsize',14)
xlabel('time', 'fontsize',14)
grid on;
zoom on;
figure(h2);
pause;

fconv=conv(f1,f2)*delta_t;
tconv=linspace(2*t_inf, 2*t_sup, 2*N-1);

h3=figure(3);
set(h3,'Position',[100 100 1200, 700])
plot(tconv,fconv,'b','linewidth',2);
title('convolution result', 'fontsize',14)
xlabel('time', 'fontsize',14)
grid on;
zoom on;
figure(h3);
pause;

% input, output and impulse response displayed together
h10=figure(10);
set(h10,'Position',[100 100 1200, 700]);
plot(t,f1,'r','linewidth',2); hold on;
plot(t,f2,'g','linewidth',2);
plot(tconv,fconv,'b','linewidth',2);
axis([ min([t tconv]) max([t tconv]), 5/4*min([f1 f2 fconv]), 5/4*max([f1 f2 fconv])]);
title('all signals together', 'fontsize',14)
legend('input signal','impulse response','output');
xlabel('time', 'fontsize',14)
grid on;
zoom on;
figure(h10);
pause;



delta_f=1/(delta_t*9*N);
B=1/(2*delta_t);
f_asc=linspace(-B,B,9*N);
f_asc2=linspace(-B,B,9*N);



pf1=abs(fftshift(fft([zeros(1,4*N), fconv, zeros(1,3*N+1)])))*delta_t;
pf2=abs(fftshift(fft([zeros(1,4*N), f1, zeros(1,4*N)])))*delta_t;
pf3=abs(fftshift(fft([zeros(1,4*N), f2, zeros(1,4*N)])))*delta_t;

hold off;
h4=figure(4);
set(h4,'Position',[100 100 1200, 700]);
plot(f_asc,pf2,'r','linewidth',2);
hold on;
plot(f_asc,pf3,'g','linewidth',2);
plot(f_asc,pf1,'b','linewidth',2);
title('absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h4);

hold off;
h5=figure(5);
set(h5,'Position',[100 100 1200, 700]);
plot(f_asc,20*log10(pf2),'r','linewidth',2);
hold on;
plot(f_asc,20*log10(pf3),'g','linewidth',2);
plot(f_asc,20*log10(pf1),'b','linewidth',2);
title('absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
ylabel('dB', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h5);
hold off;

hold off;
h6=figure(6);
set(h6,'Position',[100 100 1200, 700]);
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf2((9*N+1)/2+1:end)),'r','linewidth',2);
hold on;
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf3((9*N+1)/2+1:end)),'g','linewidth',2);
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf1((9*N+1)/2+1:end)),'b','linewidth',2);
title('absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
ylabel('dB', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h6);
hold off;

hold off;
h7=figure(7);
set(h7,'Position',[100 100 1200, 700]);
plot(f_asc,pf2/max(pf2),'r','linewidth',2);
hold on;
plot(f_asc,pf3/max(pf3),'g','linewidth',2);
plot(f_asc,pf1/max(pf1),'b','linewidth',2);
title('normalized absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h7);

hold off;
h8=figure(8);
set(h8,'Position',[100 100 1200, 700]);
plot(f_asc,20*log10(pf2/max(pf2)),'r','linewidth',2);
hold on;
plot(f_asc,20*log10(pf3/max(pf3)),'g','linewidth',2);
plot(f_asc,20*log10(pf1/max(pf1)),'b','linewidth',2);
title('normalized absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
ylabel('dB', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h8);
hold off;

hold off;
h9=figure(9);
set(h9,'Position',[100 100 1200, 700]);
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf2((9*N+1)/2+1:end)/max(pf2)),'r','linewidth',2);
hold on;
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf3((9*N+1)/2+1:end)/max(pf3)),'g','linewidth',2);
semilogx(f_asc((9*N+1)/2+1:end),20*log10(pf1((9*N+1)/2+1:end)/max(pf1)),'b','linewidth',2);
title('normalized absolute value of spectra of signals', 'fontsize',14)
xlabel('frequency', 'fontsize',14)
ylabel('dB', 'fontsize',14)
legend('input signal','impulse response','output');
grid on;
zoom on;
figure(h9);
hold off;

