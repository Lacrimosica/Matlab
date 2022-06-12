
% Spettro della porta

close all;
clear all;

display(' ')
tau=input('time length of the support of the Heaviside Pi?  T [sec]   ');
if isempty(tau)
    tau =1;
end
display(['T = ' num2str(tau,4)])
display(' ')

T=input('time length of the L2 space interval?  T0 [sec]   ');
if isempty(T)
    T =pi-1;
end
display(['T0 = ' num2str(T,4)])
display(' ')

f_max=input('maximum represented frequency?  f_max [Hz]   ');
if isempty(f_max)
    f_max =25;
end
display(['f_max = ' num2str(f_max,4)])

N=floor(T*f_max);
display(' ')
disp(['resulting number of shown components:   ', num2str(N)] )

f=linspace(-N/T,N/T,2*N+1);
f_thick=linspace(-N/T,N/T,2^8*N+1);
S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
S_thick=sqrt(T)*sin(pi*tau*f_thick)/pi./f_thick/T;
S(N+1)=tau/sqrt(T);
S_thick(2^7*N+1)=tau/sqrt(T);
% h=figure;
% stem(f, S);
% pause;
h=figure;
stem(f, S,'o');
hold on;

title('Fourier Series Components (Amplitude Spectrum)')
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')
zoom on;
grid on;
pause;
figure(h);
plot(f_thick, S_thick, 'k--');
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')
hold off;
pause;
pause;

%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%S(N+1)=tau/sqrt(T);
h=figure;
stem(f, abs(S).^2);hold on;
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')
title('Energy Spectrum')
zoom on;
grid on;
pause;
figure(h);
plot(f_thick, abs(S_thick).^2, 'k--');
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')
hold off;
pause;

display(' ')
disp(['energy of the signal:   ', num2str(tau)]);
disp(['energy carried by the shown components:   ', num2str(sum(abs(S).^2))]);


%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%S(N+1)=tau/sqrt(T);
h=figure;

t=[-512:512]/1024*T;
signal=zeros(size(t));
offset=(length(S)-1)/2;
figure;
xlabel('time (s)')
title('reconstructed signal')
plot(t, HPi2(tau,t));hold on;
axis([min(t),max(t), -0.2,1.2]);
pause;
signal=signal+1/sqrt(T)*S(offset+1);    % DC term added here
plot(t, real(signal));hold on;
pause;
hold off;
for n=1:offset
    signal=signal+1/sqrt(T)*S(n+offset+1)*exp(j*n*2*pi*t/T)+1/sqrt(T)*S(offset-n+1)*exp(-j*n*2*pi*t/T);
    plot(t, HPi2(tau,t));hold on;
    % plot(t, real(signal));hold off;
    plot(t, real(signal));hold on;
    axis([min(t),max(t), -0.2,1.2]);
    pause;
end;
pause;
hold off;
xlabel('time (s)')
title('reconstructed signal')
plot(t, HPi2(tau,t));hold on;
axis([min(t),max(t), -0.2,1.2]);
pause;
plot(t, real(signal));hold off;
zoom on;
grid on;


display(' ')
disp(['energy of the signal:   ', num2str(tau)]);
disp(['energy carried by the shown components:   ', num2str(sum(abs(S).^2))]);

