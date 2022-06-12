
% Spettro della porta

close all;
clear all;

display(' ')
T=input('time length of the support of the Heaviside Pi?  T [sec]   ');
if isempty(T)
    T =1;
end
display(['T = ' num2str(T,4)])
display(' ')

T0=input('time length of the L2 space interval?  T0 [sec]   ');
if isempty(T0)
    T0 =pi-1;
end
display(['T0 = ' num2str(T0,4)])
display(' ')

f_max=input('maximum represented frequency?  f_max [Hz]   ');
if isempty(f_max)
    f_max =10;
end
display(['f_max = ' num2str(f_max,4)])

N=floor(T0*f_max);
display(' ')
disp(['resulting number of shown components:   ', num2str(N)] )





f=linspace(-N/T0,N/T0,2*N+1);
f_thick=linspace(-N/T0,N/T0,2^8*N+1);
S=sin(pi*T*f)/pi./f;
S_thick=sin(pi*T*f_thick)/pi./f_thick;
S(N+1)=T;
S_thick(2^7*N+1)=T;



t=[-512:512]/1024*T0;
signal=zeros(size(t));
offset=(length(S)-1)/2;
figure;
xlabel('time (s)')
title('reconstructed signal')
plot(t, HPi2(T,t));hold on;
axis([min(t),max(t), 1.4*min(S),1.2*max(S)]);
pause;
% signal=signal+1/sqrt(T)*S(offset+1);    % DC term added here
% plot(t, real(signal));
% pause;



% h=figure;
% stem(f, S);
% pause;
h=figure;
plot(f, S,'o','markersize',4);
gcah=gca;
gcah.FontSize=12;
axis([min(f),max(f), 1.4*min(S),1.2*max(S)]);
%stem(f, S,'o');
hold on;

title('Fourier Basis Projections s_{n}')
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')
zoom on;
grid on;
axis([])
pause;
figure(h);

plot(f_thick, S_thick, 'k--', 'LineWidth', 0.2);
xlabel('frequency fn = n/T_{0} = (n f_{0})   [Hz]')

hold off;

return

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
disp(['energy of the signal:   ', num2str(T)]);
disp(['energy carried by the shown components:   ', num2str(sum(abs(S).^2))]);


%S=sqrt(T)*sin(pi*tau*f)/pi./f/T;
%S(N+1)=tau/sqrt(T);
h=figure;

t=[-512:512]/1024*T0;
signal=zeros(size(t));
offset=(length(S)-1)/2;
figure;
xlabel('time (s)')
title('reconstructed signal')
plot(t, HPi2(T,t));hold on;
axis([min(t),max(t), -0.2,1.2]);
pause;
signal=signal+1/sqrt(T0)*S(offset+1);    % DC term added here
plot(t, real(signal));hold on;
pause;
hold off;
for n=1:offset
    signal=signal+1/sqrt(T0)*S(n+offset+1)*exp(j*n*2*pi*t/T0)+1/sqrt(T0)*S(offset-n+1)*exp(-j*n*2*pi*t/T0);
    plot(t, HPi2(T,t));hold on;
    % plot(t, real(signal));hold off;
    plot(t, real(signal));hold on;
    axis([min(t),max(t), -0.2,1.2]);
    pause;
end;
pause;
hold off;
xlabel('time (s)')
title('reconstructed signal')
plot(t, HPi2(T,t));hold on;
axis([min(t),max(t), -0.2,1.2]);
pause;
plot(t, real(signal));hold off;
zoom on;
grid on;


display(' ')
disp(['energy of the signal:   ', num2str(T)]);
disp(['energy carried by the shown components:   ', num2str(sum(abs(S).^2))]);

