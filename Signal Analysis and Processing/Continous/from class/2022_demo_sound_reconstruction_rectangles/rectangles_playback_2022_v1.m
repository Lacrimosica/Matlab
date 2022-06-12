

clear all
close all
clc
load results2_18.mat

display(['time-width of rectangles: ', num2str(t0*1e3), '  ms']);
display(['signal energy: ', num2str(signal_energy)]);
display(['approximant energy: ', num2str(approximant_energy)]);
error_energy=signal_energy-approximant_energy;
display(['approximant-to-signal ratio: ', num2str(approximant_energy/signal_energy*100), '  %']);
display(['noise-to-signal ratio: ', num2str(error_energy/signal_energy*100), '  %']);
display(['signal-to-noise ratio dB: ', num2str(10*log10(signal_energy/error_energy))]);

fs=2^20;   % 1048576 kHz
f0=1/t0;
Q=round(fs/f0);
ts=t0/Q;
t=[0:ts:t_end-ts];

original_signal=test_signal2(t+ts/2);
approximate_signal=s(floor(t/t0)+1)/sqrt_t0;

figure; plot(t,original_signal);
figure; plot(t,original_signal, t,approximate_signal);
xlabel('time (s)', 'fontsize', 16);
title(['error signal vs original signal energy: ' num2str((signal_energy-approximant_energy)/signal_energy*100), '  %'], 'fontsize', 16);
gcah=gca;
gcah.FontSize=16;


% sum(original_signal.^2)*ts
% sum(approximate_signal.^2)*ts
% (sum(original_signal.^2)-sum(approximate_signal.^2))/sum(original_signal.^2)


pause;

fs=2^17;   % 131072 kHz
f0=1/t0;
Q=(fs/f0);
ts=t0/Q;
t=[0:ts:t_end-ts];

original_signal=test_signal2(t+ts/2);
approximate_signal=s(floor(t/t0)+1)/sqrt_t0;

sound(original_signal,fs);
pause;
sound(approximate_signal,fs);
pause;
sound(original_signal-approximate_signal',fs);

