
% unilateral step function

close all;
clear all;

display(' ')
a_const=input('exponential decay coefficient a (s^-1) ?  [default 1/2] ');
if isempty(a_const)
    a_const =1/2;
end
display(['the exponential decay coefficient a is:  ' num2str(a_const,3) '  (s^-1)'])

display(' ')
T0=input('time length of the interval T0 (s) ?  [default 10]   ');
if isempty(T0)
    T0 =10;
end
f0= 1/T0;
display(['the time length of the interval T0 is:  ' num2str(T0,3) '  (s)'])

display(' ')
display(['the interval is [' num2str(-T0/2,3) ' , ' num2str(T0/2,3) ' ]'])
display(['the frequency step f0 is:  ' num2str(f0,3)])

display(' ')
f_max=input('maximum plotted frequency  f_max (Hz) ?  [default 10]  ');
if isempty(f_max)
    f_max =10;
end
N=floor( f_max / f0 );
display(['the maximum plotted frequency is:  ' num2str(f_max,3) '  (Hz)'])
display(['the resulting number of shown spectral components is:  ', num2str(N)] )

display(' ')
f_cut_off=input('maximum frequency used for reconstruction  f_cut_off (Hz) ?  [default 20]   ');
if isempty(f_cut_off)
    f_cut_off =5;
end
display(['the maximum frequency used for reconstruction is:  ', num2str(f_cut_off,3) '  (Hz)'] )


% thickening factor, or representing time-continuous 
% or frequency-continuous plots
NFF=2^8;

n=-N:N;
nf0=linspace(-N*f0,N*f0,2*N+1);
t=linspace(-T0/2,T0/2,2*N+1);
%t=linspace(0,T,2*N+1);
f_thick=linspace(-N*f0,N*f0,NFF*2*N+1);
t_thick=linspace(-T0/2,T0/2,NFF*2*N+1);
%t_thick=linspace(0,T,NFF*2*N+1);
%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%
% time-continuous signal
unilat_exp=exp(-a_const*t_thick).*u(t_thick);
% Fourier components
sn=1/sqrt(T0)*(1-exp(-a_const*T0/2).*((-1).^n))./(a_const+1i*2*pi*nf0);
% frequency-continuous spectrum
S_thick=1/sqrt(T0)*(1-exp(-a_const*T0/2).*exp(-1i*pi*f_thick/f0))./(a_const+1i*2*pi*f_thick);

% h=figure;
% stem(f, S);
% pause;

ht=figure;
plot(t_thick,unilat_exp,'linewidth',2);
axis([min(t_thick), max(t_thick), ...
    (max(unilat_exp)+min(unilat_exp))/2-abs(max(unilat_exp)-min(unilat_exp))*1.25/2,...
    (max(unilat_exp)+min(unilat_exp))/2+abs(max(unilat_exp)-min(unilat_exp))*1.25/2]);
title('unilateral exponential signal')
xlabel('time  [s]')
zoom on;
grid on;
figure(ht);
pause;

% plotting spectra
hs=figure;
stem(nf0, real(sn),'*g','markersize',8);
hold on;
title('Fourier Series Components (Amplitude Spectrum)')
xlabel('frequency f_n = n/T_{0} = (n f0)   [Hz]')
zoom on;
grid on;
legend(' real part')
pause;
stem(nf0, imag(sn),'or','markersize',8);
legend(' real part', 'imaginary part' )
pause;
figure(hs);
plot(f_thick, real(S_thick), 'g--');
pause;
plot(f_thick, imag(S_thick), 'r--');
hold off;
component_axis=axis;
pause;

%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%S(N+1)=tau/sqrt(T);
hss=figure;
stem(nf0, abs(sn).^2,'db');hold on;
xlabel('frequency fn = n/T_{0} = (n f0)   [Hz]')
title('Energy Spectrum')
zoom on;
grid on;
pause;
figure(hss);
plot(f_thick, abs(S_thick).^2, 'b--');
xlabel('frequency fn = n/T_{0} = (n f0)  [Hz]')
hold off;
pause;

% raised_thick=root_raised_cosine(t_thick,alpha,tau,0,0);
% disp(['energy of the signal:   ', num2str(sum(abs(raised_thick).^2)*T/(2^12*N+1))]);
% disp(['energy carried by the shown components:   ', num2str(sum(abs(sn).^2))]);

%%% reconstruction
% apply frequency limitation
N_cut_off=floor(T0*f_cut_off);
% indices used in the summation
n_app=-N_cut_off:N_cut_off;
disp([' resulting number of reconstruction components (positive):   ', num2str(N_cut_off)] )
hs_app=figure;
nf0=linspace(-N_cut_off*f0,N_cut_off*f0,2*N_cut_off+1);
sn_app=1/sqrt(T0)*(1-exp(-a_const*T0/2).*((-1).^n_app))./(a_const+1i*2*pi*nf0);
stem(nf0, real(sn_app),'*g','markersize',8);hold on;
stem(nf0, imag(sn_app),'or'),'markersize',8;
axis(component_axis);
hold off;
title('frequency components used for reconstruction')
pause;

% actually reconstructing
unilat_exp_app=zeros(1,length(t_thick));
for k=1:2*N_cut_off+1
    unilat_exp_app=unilat_exp_app+real(sn_app(k).*exp(1i*2*pi*t_thick/T0*(k-N_cut_off-1))/sqrt(T0));
end;

figure(ht);hold on;pause;
plot(t_thick,unilat_exp,'b','linewidth',2);
plot(t_thick,unilat_exp_app,'r','linewidth',2);
title('unilateral exponential signal (blue) and approximation (red)')
xlabel('time  (s)')
axis([min(t_thick), max(t_thick), ...
    (max(unilat_exp_app)+min(unilat_exp_app))/2-abs(max(unilat_exp_app)-min(unilat_exp_app))*1.25/2,...
    (max(unilat_exp_app)+min(unilat_exp_app))/2+abs(max(unilat_exp_app)-min(unilat_exp_app))*1.25/2]);
zoom on;
grid on;

% energy calculations
% energy of original signal (from Parseval's rule, be better if analytical)
Es=(1-exp(-2*a_const*T0/2))/(2*a_const);
% energy of reconstructed signal from Parseval's rule
E_comp_app=sum(abs(sn_app).^2)
% error energy in percentage of signal energy
E_ratio_perc=(Es-E_comp_app)/Es*100

display(' ')
display(['the energy of the signal is: ' num2str(Es,5) ])
display(['the energy of approximant is: ' num2str(E_comp_app,5) ])
display(['the error energy in percentage of the signal energy is: ' num2str(E_ratio_perc,5) ])


