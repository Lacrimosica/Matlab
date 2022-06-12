%% Fourier Transforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
% setting the time interval for displaying the time-domain signal
t=[-25:0.02:25];
% setting the frequency points where the Fourier transform is calculated
f=[-10:0.01:10];
% counting how many time or frequency points there are
Nt=numel(t);
Nf=numel(f);
% string initialization for print out of calculation progress
strPerc = 0; strMsg = sprintf('\n percentage executed %d/%d',strPerc);
%
% ideal lowpass filter
B=1; H= @(f) HPi(2*B,f);
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H(f)),'r',f,imag(H(f)),'b','linewidth',2);
Dr=1.25*(max(max(real(H(f))),max(imag(H(f))))-min(min(real(H(f))),min(imag(H(f)))));
mr=(max(max(real(H(f))),max(imag(H(f))))+min(min(real(H(f))),min(imag(H(f)))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the ideal lowpass filter transfer function in frequency domain');
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% impulse response of the ideal lowpass filter:
% ideal lowpass filter
B=1; h= @(t) 2*B*sinc(2*B*t);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');
plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title('plot of the impulse response of the ideal lowpass filter in time domain'); grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% truncating the ideal lowpass filter impulse response
% set here the total duration of the truncated filter impulse response
trunc_filt_duration=10/2/B;
h= @(t) 2*B*sinc(2*B*t).*HPi(trunc_filt_duration,t);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title(['plot of the of the ideal lowpass filter IMPULSE RESPONSE TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']); grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% calculating the transfer function of the truncated ideal lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -10, 10);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the ideal lowpass filter TRANSFER FUNCTION when the impulse response is TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']); grid on;
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% truncating the ideal lowpass filter impulse response WIDER
% set here the total duration of the truncated filter impulse response
trunc_filt_duration=30/2/B;
h= @(t) 2*B*sinc(2*B*t).*HPi(trunc_filt_duration,t);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title(['plot of the of the ideal lowpass filter IMPULSE RESPONSE TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% calculating the transfer function of the truncated ideal lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -20, 20);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the ideal lowpass filter TRANSFER FUNCTION when the impulse response is TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%

%% looking at raised-cosine transfer functions
B=1; alpha=0.2; H= @(f) rho(f,alpha,2*B,0,0);
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H(f)),'r',f,imag(H(f)),'b','linewidth',2);
Dr=1.25*(max(max(real(H(f))),max(imag(H(f))))-min(min(real(H(f))),min(imag(H(f)))));
mr=(max(max(real(H(f))),max(imag(H(f))))+min(min(real(H(f))),min(imag(H(f)))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the raised-cosine lowpass filter transfer function with roll-off ' num2str(alpha)]);
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
B=1; alpha=0.2; h= @(t) 2*B*sinc(2*B*t).*cos(2*B*pi*alpha*t)./(1-4*alpha^2*t.^2*4*B^2);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title(['plot of the impulse response of the raised-cosine lowpass filter with roll-off ' num2str(alpha)]);
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% truncating the raised-cosine lowpass filter impulse response
% set here the total duration of the truncated filter impulse response
trunc_filt_duration=10/2/B;
B=1; alpha=0.2; h= @(t) 2*B*sinc(2*B*t).*cos(2*B*pi*alpha*t)./(1-4*alpha^2*t.^2*4*B^2).*HPi(trunc_filt_duration,t);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title(['plot of the of the raised-cosine lowpass filter IMPULSE RESPONSE TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% calculating the transfer function of the truncated raised-cosine lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -10, 10);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the raised-cosine lowpass filter TRANSFER FUNCTION when the impulse response is TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% truncating the raised-cosine lowpass filter impulse response to a WIDER time
% set here the total duration of the truncated filter impulse response
trunc_filt_duration=30/2/B;
B=1; alpha=0.2; h= @(t) 2*B*sinc(2*B*t).*cos(2*B*pi*alpha*t)./(1-4*alpha^2*t.^2*4*B^2).*HPi(trunc_filt_duration,t);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title(['plot of the of the raised-cosine lowpass filter IMPULSE RESPONSE TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']); 
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% calculating the transfer function of the longer truncated raised-cosine lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -10, 10);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the raised-cosine lowpass filter TRANSFER FUNCTION when the impulse response is TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% plot of the absolute value of the Fourier Transform
figure('windowstate','maximized');plot(f,abs(H),'k','linewidth',2);
Dr=1.25*(max(abs(H))-min(abs(H)));
mr=(max(abs(H))+min(abs(H)))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the ABSOLUTE VALUE of the raised-cosine lowpass filter TRANSFER FUNCTION when the impulse response is TRUNCATED to ' num2str(trunc_filt_duration) ' (s)']);
xlabel('frequency [Hz]');legend('absolute value');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
%% making causal the truncated impulse response of the raised-cosine filter 
B=1; alpha=0.2; h= @(t) 2*B*sinc(2*B*(t-trunc_filt_duration/2)).*cos(2*B*pi*alpha*(t-trunc_filt_duration/2))...
    ./(1-4*alpha^2*(t-trunc_filt_duration/2).^2*4*B^2).*HPi(trunc_filt_duration,t-trunc_filt_duration/2);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title('plot of the impulse response in time domain'); 
title(['plot of the raised-cosine lowpass filter impulse response truncated to ' num2str(trunc_filt_duration) ' (s) and MADE CAUSAL']); 
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% calculating the transfer function of the longer truncated raised-cosine lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -5, 20);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the raised-cosine filter transfer function with truncated and made-causal impulse response of duration ' num2str(trunc_filt_duration) ' (s)']); 
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
pause;
%
% plot of the absolute value of the Fourier Transform
figure('windowstate','maximized');plot(f,abs(H),'k','linewidth',2);
Dr=1.25*(max(abs(H))-min(abs(H)));
mr=(max(abs(H))+min(abs(H)))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title(['plot of the ABSOLUTE VALUE of the raised-cosine filter transfer function with truncated and made-causal impulse response of duration ' num2str(trunc_filt_duration) ' (s)']); 
xlabel('frequency [Hz]');legend('absolute value');
grid on;
current_gca=gca;current_gca.FontSize=16;
zoom on;
%
return

%% THIS PART NOT USED THIS YEAR
%% "Sample and Hold" impulse response
%
B=1; h= @(t) HPi(1/2/B,t-1/4/B);
% plot of the real and imaginary part of the impulse response in time-domain
figure('windowstate','maximized');plot(t,real(h(t)),t,imag(h(t)),'linewidth',2);
Dr=1.25*(max(max(real(h(t))),max(imag(h(t))))-min(min(real(h(t))),min(imag(h(t)))));
mr=(max(max(real(h(t))),max(imag(h(t))))+min(min(real(h(t))),min(imag(h(t)))))/2;
axis([t(1),t(end),mr-Dr/2,mr+Dr/2]);
xlabel('time [s]');legend('real part', 'imaginary part');
title('plot of the impulse response in time domain'); grid on;
pause;
%
% calculating the transfer function of the longer truncated raised-cosine lowpass filter
% calculation loop
clear H;
for n=1:Nf;
    % 
    %%% actual Fourier Transform calculation takes place here
    H(n)=integral(@(t) h(t) .* exp(-1i*2*pi*f(n)*t) , -5, 20);
    %%%
    %
    % displaying progress on console
    strErase = sprintf(repmat('\b',1,numel(strMsg)));
    strMsg = sprintf('\n percentage executed %i',round(n/Nf*100));
    fprintf([strErase strMsg]);
    %
end
%
% plot of the real and imaginary part of the Fourier Transform
figure('windowstate','maximized');plot(f,real(H),'r',f,imag(H),'b','linewidth',2);
Dr=1.25*(max(max(real(H)),max(imag(H)))-min(min(real(H)),min(imag(H))));
mr=(max(max(real(H)),max(imag(H)))+min(min(real(H)),min(imag(H))))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the transfer function in frequency domain');
xlabel('frequency [Hz]');legend('real part', 'imaginary part');
grid on; 
pause;
%
% plot of the absolute value of the Fourier Transform
figure('windowstate','maximized');plot(f,abs(H),'k','linewidth',2);
Dr=1.25*(max(abs(H))-min(abs(H)));
mr=(max(abs(H))+min(abs(H)))/2;
axis([f(1),f(end),mr-Dr/2,mr+Dr/2]);
title('plot of the absolute value of the transfer function in frequency domain');
xlabel('frequency [Hz]');legend('absolute value');
grid on;

