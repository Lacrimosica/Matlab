%% Fourier Transforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% string initialization for print out of calculation progress
strPerc = 0; strMsg = sprintf('\n percentage executed %d/%d',strPerc);
% switch inizialization (see below for usage)
elem_sig_switch=['off'];
%
%these initializations are meant to inhibit graphic display for periodic
%signal, when they are not used
nmax=-1;nmin=0;T0=0;
%
%% signal definitions (uncomment the desired one, one signal per line)
% T=1/2; s= @(t) HPi(T,t);
% T=2; a=0.25; td=2; s= @(t) HPi(T,a*(t-td));
% T=1; s= @(t) HPi(T,t+T/2)-HPi(T,t-T/2);
a=1; s= @(t) u(t).*exp(-a*t) - u(-t).*exp(a*t);
% T=1; td=0.2; s= @(t) HPi(T,t-td);
% T=1; a=2; s= @(t) HPi(T,a*t);
% T=1; td=3; s= @(t) HLambda(T,t-td);
% a=1; s= @(t) u(t).*exp(-a*t);
% a=1; s= @(t) u(-t).*exp(a*t);
% a=1/2; s= @(t) u(t).*t.*exp(-a*t);
% a=1; s= @(t) u(t).*exp(-a*t)+u(-t).*exp(a*t);
% a=1; s= @(t) u(t).*exp(-a*t)-u(-t).*exp(a*t);
% a=1; n=1; s= @(t) u(t).*(t.^n).*exp(-a*t);
% a=1; s= @(t) exp(-a*abs(t));
% T=1; alpha=0.5; s= @(t) rho(t,alpha,T,0,0);
% alpha=10; a=1; s= @(t) u(alpha*t).*(alpha*t).*exp(-a*alpha*t);
% alpha=10; a=1; n=2; s= @(t) u(alpha*t).*(alpha*t.^n).*exp(-a*alpha*t);
% T=1; td=1; s= @(t) HPi(T,t-td);
% T=1; td=1; s= @(t) t.*HPi(T,t-td);
% T=1; s= @(t) t.*HLambda(T,t);
% a=1; alpha=4; s= @(t) 1./(a+1i*2*pi*t*alpha);
% T0=1; s= @(t) (cos(2*pi*t/T0)+1).*HPi(T0,t)/T0;
% T=1; f0=5; s= @(t) exp(j*2*pi*t*f0).*HPi(T,t);
% T=1; f0=10; s= @(t) cos(2*pi*f0*t).*HPi(T,t)/T;
% s= @(t) HLambda(1,t-1/2)+HLambda(1,t+1/2)+1/2*HLambda(1,2*t);
%
%% periodic signals (uncomment the desired one, one signal per line, and uncomment the "periodicizer")
% define the elementary signal "q(t)"
% T=1/2; q = @(t) HPi(T,t);
% T=1/4; q = @(t) HPi(T,t-T/2);
% T=1/2; q = @(t) HPi(T/2,t-T/4)+HPi(T/2,t-7/4*T);
% T=1/2; q= @(t) HPi(T,t+T/2)-HPi(T,t-T/2);
% sigma=1/3; q = @(t) exp((-t.^2)/2/sigma^2);
% a=1; q = @(t) exp(-a*t).*u(t);
% T=1/10; alpha=0.3; q= @(t) rho(t,alpha,T,0,0);
% T0=1; T=1/2; q = @(t) 2*HPi(T,t)-HPi(T0,t);
% T=1/17; q = @(t) HPi(T,t);
% T=1/2; q = @(t) HPi(T,t+T/2)-HPi(T,t-T/2);
%
%% define a (quasi)-periodic version of it, where T0 is the period, 
% and nmin, nmax are the summation index limits 
% T0=1; nmin=-10; nmax=10; s= @(t) periodicizer(q,T0,nmin,nmax,t);
%% set the switch to "on" if you want to see the individual elementary
% signals making up the periodic signal
elem_sig_switch=['on'];
%
% setting the time interval for calculating the Fourier Transform
% for accurate calculation the signal should be zero outside this interval
t0=-11; t1=11;
%
% Setting the time interval for displaying the signal
t=[-15:0.01:15];
% setting the frequency points where the Fourier transform is calculated
f=[-20:0.02:20];
% counting how many time or frequency points there are
Nt=numel(t);
Nf=numel(f);
%
%% calculation loop
for n=1:Nf;
    % 
    % actual calculation takes place here
    S(n)=integral(@(t) s(t) .* exp(-1i*2*pi*f(n)*t) , t0, t1, 'RelTol',1e-4);
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
%
%% outputs and figures
disp(' ');
%
figure('windowstate','maximized');
%

% plotting the individual signals making up the periodic signal
if strcmp(elem_sig_switch,['on'])
    for n=1:(nmax-nmin+1)
        plot(t,real(q(t-(n+nmin-1)*T0)),'--','linewidth',1);hold on;
        plot(t,imag(q(t-(n+nmin-1)*T0)),'-.','linewidth',1);
    end;
    Dr=1.25*(max(max(real(s(t))),max(imag(s(t))))-min(min(real(s(t))),min(imag(s(t)))));
    mr=(max(max(real(s(t))),max(imag(s(t))))+min(min(real(s(t))),min(imag(s(t)))))/2;
    axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);

    title('plot of the signal in time domain'); grid on;
    current_gca=gca;current_gca.FontSize=16;
    zoom on;
    pause;
end;
% plotting the overall signal
plot(t,real(s(t)),t,imag(s(t)),'linewidth',2);
Dr=1.25*(max(max(real(s(t))),max(imag(s(t))))-min(min(real(s(t))),min(imag(s(t)))));
mr=(max(max(real(s(t))),max(imag(s(t))))+min(min(real(s(t))),min(imag(s(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
title('plot of the signal in time domain'); grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
xlabel('time [s]');legend('real part', 'imaginary part');
pause;
%

figure('windowstate','maximized');plot(f,real(S),f,imag(S),'linewidth',2);
Dr=1.25*(max(max(real(S)),max(imag(S)))-min(min(real(S)),min(imag(S))));
mr=(max(max(real(S)),max(imag(S)))+min(min(real(S)),min(imag(S))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the signal in frequency domain');
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
figure('windowstate','maximized');plot(f,abs(S),'k','linewidth',2);
Dr=1.25*(max(abs(S))-min(abs(S)));
mr=(max(abs(S))+min(abs(S)))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the absolute value of the signal in frequency domain');
xlabel('frequency [Hz]');legend('absolute value');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
figure('windowstate','maximized');plot(f,abs(S).^2,'k','linewidth',2);
Dr=1.25*(max(abs(S).^2)-min(abs(S).^2));
mr=(max(abs(S).^2)+min(abs(S).^2))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the energy spectrum of the signal in frequency domain');
xlabel('frequency [Hz]');legend('energy spectrum');
grid on;
zoom on;
current_gca=gca;current_gca.FontSize=16;

pause;
close all;

