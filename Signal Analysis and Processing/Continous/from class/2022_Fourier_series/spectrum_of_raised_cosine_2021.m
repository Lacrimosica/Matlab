
% Raised Cosine

close all;
clear all;

tau=input('duration at half height?  T0 [sec]   ');
if isempty(tau)
    tau =1
end
alpha=input('roll-off factor?  alpha   ');
if isempty(alpha)
    alpha =0.5
end
T=input('time length of the L2 space interval?  T [sec]   ');
if isempty(T)
    T =pi*tau
end
f_max=input('maximum represented frequency?  f_max [Hz]   ');
if isempty(f_max)
    f_max =25
end
f_cut_off=input('maximum frequency used for reconstruction ?  f_cut_off [Hz]   ');
if isempty(f_cut_off)
    f_cut_off =10
end

N=floor(T*f_max);
disp([' resulting number of shown components:   ', num2str(N)] )

f=linspace(-N/T,N/T,2*N+1);
t=linspace(-T/2,T/2,2*N+1);
f_thick=linspace(-N/T,N/T,2^10*N+1);
t_thick=linspace(-T/2,T/2,2^10*N+1);
%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
S=sqrt(T)*sin(pi*tau*f)/pi./f/T .* cos(pi*alpha*f*tau)./(1-4*alpha^2.*f.^2*tau^2);
%S_thick=sqrt(T)*sin(pi*tau*f_thick)/pi./f_thick/T;
S_thick=sqrt(T)*sin(pi*tau*f_thick)/pi./f_thick/T .* cos(pi*alpha*f_thick*tau)./(1-4*alpha^2.*f_thick.^2*tau^2);
% explicitly calculating the DC component
S(N+1)=tau/sqrt(T);
S_thick(2^9*N+1)=tau/sqrt(T);
raised=root_raised_cosine(t_thick,alpha,tau,0,0);
% h=figure;
% stem(f, S);
% pause;

ht=figure;
plot(t_thick,raised,'linewidth',2);
axis([min(t_thick), max(t_thick), ...
    (max(raised)+min(raised))/2-abs(max(raised)-min(raised))*1.25/2,...
    (max(raised)+min(raised))/2+abs(max(raised)-min(raised))*1.25/2]);
title('raised cosine signal')
xlabel('time  [s]')
zoom on;
grid on;
figure(ht);
pause;


hs=figure;
stem(f, S,'o');
hold on;
title('Fourier Series Components (Amplitude Spectrum)')
xlabel('frequency fn = n/T = (n f0)   [Hz]')
zoom on;
grid on;
pause;
figure(hs);
plot(f_thick, S_thick, 'k--');
xlabel('frequency fn = n/T = (n f0)   [Hz]')
hold off;
pause;

%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%S(N+1)=tau/sqrt(T);
hss=figure;
stem(f, abs(S).^2);hold on;
xlabel('frequency fn = n/T = (n f0)   [Hz]')
title('Energy Spectrum')
zoom on;
grid on;
pause;
figure(hss);
plot(f_thick, abs(S_thick).^2, 'k--');
xlabel('frequency fn = n/T = (n f0)  [Hz]')
hold off;
pause;

raised_thick=root_raised_cosine(t_thick,alpha,tau,0,0);
disp(['energy of the signal:   ', num2str(sum(abs(raised_thick).^2)*T/(2^12*N+1))]);
disp(['energy carried by the shown components:   ', num2str(sum(abs(S).^2))]);

%%% reconstruction
% apply frequency limitation
N_cut_off=floor(T*f_cut_off);
disp([' resulting number of reconstruction components:   ', num2str(N_cut_off)] )
figure(hs); hold on;
f=linspace(-N_cut_off/T,N_cut_off/T,2*N_cut_off+1);
S_app=sqrt(T)*sin(pi*tau*f)/pi./f/T .* cos(pi*alpha*f*tau)./(1-4*alpha^2.*f.^2*tau^2);
S_app(N_cut_off+1)=tau/sqrt(T);
stem(f, S_app,'or');
hold off;
title('red: frequency components used for reconstruction')
pause;
% actually reconstructing
raised_app=zeros(1,length(t_thick));
for k=1:2*N_cut_off+1
    raised_app=raised_app+real(S_app(k).*exp(1i*2*pi*t_thick/T*(k-N_cut_off-1))/sqrt(T));
end;

figure(ht);hold on;pause;
plot(t_thick,raised_app,'r','linewidth',2);
title('raised cosine signal (blue) and approximation (red)')
xlabel('time  [s]')
zoom on;
grid on;
    



