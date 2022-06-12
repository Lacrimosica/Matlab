clear all;
close all;

load vasca x
% load classical x

% J.S.Bach; Partita E major, Gavotte en rondeau (excerpt)
% Sirkka Väisänen, violin. 

M=330750;
Fs=44100;
%y=x(1:M,1)+x(1:M,2);
%ym=mean(y);
%y=y-ym;

% zero padding by two seconds at beginning and 2 seconds at end of song
Nds=0;  % MUST BE AN INTEGER
y=zeros(Nds*Fs+M,1);
% stereo for "classical"
% y(Nds/2*Fs+1:M+Nds/2*Fs)=x(1:M,1)+x(1:M,2);
% mono for vasca
y(Nds/2*Fs+1:M+Nds/2*Fs)=x(1:M);
% the new number of time samples is
M=M+Nds*Fs;

clear x;

% volume control
y=y*0.5;
% killing DC value
y=y-mean(y);

% playing song
pause;
sound(y,Fs);
%pause;
t0=1/Fs;
t=[0:1:M-1]*t0;
T=M*t0;
f0=Fs/M;
display(' ')
display(['The total duration of the signal T is :' num2str(T)]);
display(['The fundamental frequency f0 is :' num2str(f0)]);
display(['The number of considered positive-frequency components N is :' num2str(M/2)]);
display(['The maximum represented positive-frequency (N*f0) is :' num2str(M/2*f0)]);
display(' ')

% plotting the signal

%% before plotting the signal, it is oversampled
% oversampling factor
Nos=8;
yf=fft(y);
yfover=[yf(1:M/2) ; (0+1i*zeros(1,(Nos-1)*M))' ; yf(M/2+1:M)];
yover=Nos*ifft(yfover);
tover=[0:1:Nos*M-1]*t0/Nos;

figure('windowstate','maximized');
%plot(t,y);
title('SIGNAL PLOT');
xlabel('time [s]');
zoom on; 
%hold on;
plot(tover,real(yover));
pause;
%%


f=[-M/2:1:(M/2-1)]*(Fs/2)/(M/2);
% calcolo dello spettro del segnale
yf=fft(y);
yfdisplay=fftshift(yf);
% plot dello spettro del segnale
% plot(f,abs(yfdisplay));
figure('windowstate','maximized');
stem(f,abs(yfdisplay),'marker','none');
title('SIGNAL SPECTRUM');
xlabel('frequency [Hz]');
ylabel('absolute value of Fourier components');
grid on;
pause;




figure('windowstate','maximized');
plot(f,20*log10(abs(yfdisplay)));
% stem(f,20*log10(abs(yfdisplay)),'marker','none');
%semilogy(f,abs(yfdisplay));
title('SIGNAL SPECTRUM');
xlabel('frequency [Hz]');
ylabel('absolute value of Fourier components, [dB]');
grid on;

% evaluating the energy of the signal
% in frequency:
Ef=sum(yf.*conj(yf))/M/Fs;
% in time:
Et=sum(y.*conj(y))/Fs;

display(['The signal energy is (time integral): ' num2str(Et)]) 
display(['The signal energy is (sum of components abs-squared): ' num2str(Ef)]) 
display(' ')

pause;

keep_on_going=1;

while keep_on_going

% killing the high end of the band
highcut=[];
while isempty(highcut)
highcut=input('input cut-off frequency:  ');
display(' ')
end;
%highcut=15000;
maxfsample=round(highcut/f0);
if maxfsample>M/2
    maxfsample=M/2;
end;
N=M/2-maxfsample;
yf1=[yf(1:M/2-N);zeros(2*N+1,1)+j*zeros(2*N+1,1);yf(M/2+N+2:M)];

hold on;
yf1display=fftshift(yf1);
figure('windowstate','maximized');
% plot(f,abs(yfdisplay));
stem(f,abs(yfdisplay),'marker','none');
hold on;
% plot(f,abs(yf1display),'r');
stem(f,abs(yf1display),'r','marker','none');
title('SIGNAL SPECTRUM WITH CUT-OFF');
xlabel('frequency [Hz]');
ylabel('absolute value of Fourier components');
zoom on;
grid on;
pause;
hold off;

% figure;
% plot(f,abs(yfdisplay).^2);
% grid on;
% hold on;
% pause;
% plot(f,abs(yf1display).^2,'r');
% zoom on;
% pause;
% hold off;

figure('windowstate','maximized');
plot(f,20*log10(abs(yfdisplay)));
hold on;
plot(f,20*log10(abs(yf1display)),'r');
title('SIGNAL SPECTRUM WITH CUT-OFF');
xlabel('frequency [Hz]');
ylabel('absolute value of Fourier components, [dB]');
grid on;
zoom on;
pause;
hold off;

% frequency band-limited signal
y1=ifft(yf1);

% figure('windowstate','maximized');
% plot(t,y);
% hold on;
% plot(t,real(y1),'r');
% title('SIGNAL PLOT AFTER FREQUENCY CUT-OFF');
% xlabel('time [s]');
% zoom on;
% pause;
% hold off;

%% before plotting the signal, it is oversampled
% oversampling factor
Nos=8;
yf=fft(y);
yfover=[yf(1:M/2) ; (0+1i*zeros(1,(Nos-1)*M))' ; yf(M/2+1:M)];
yover=Nos*ifft(yfover);
yf1=fft(y1);
y1fover=[yf1(1:M/2) ; (0+1i*zeros(1,(Nos-1)*M))' ; yf1(M/2+1:M)];
y1over=Nos*ifft(y1fover);
tover=[0:1:Nos*M-1]*t0/Nos;

figure('windowstate','maximized');
plot(tover,real(yover));
hold on;
plot(tover,real(y1over),'r');
title('SIGNAL PLOT AFTER FREQUENCY CUT-OFF');
xlabel('time [s]');
zoom on;
pause;
hold off;

%%

sound(real(y1),Fs);


% Energy of the signal:
% evaluating the energy of the signal
% in frequency:
Ef=sum(yf1.*conj(yf1))/M/Fs;
Et=sum(y1.*conj(y1))/Fs;
Eerr=sum((y-y1).*conj(y-y1))/Fs;
display(['The approximant signal energy is (time integral): ' num2str(Et)]) ;
display(['The approximant signal energy is (sum of components abs-squared): ' num2str(Ef)]) ;
display(['The energy of the error signal is (time integral): ' num2str(Eerr)]);
display(['Check of the sum of the two: ' num2str(Ef+Eerr)]);
display(['The error energy is ' num2str(Eerr/Ef*100) '% of the original signal']);
display(' ')

go_on=[];
while isempty(go_on)
go_on=input('do you want to continue? (y/n)  ','s');
display(' ')
end;
keep_on_going=~strcmp(go_on,'n');

end