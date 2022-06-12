clear all
close all
clc


s=tf('s')
H= 1/((s+2) * (s+10));
w0=0.5;

[magnitude, phase] = bode(H, w0)

%%it will output phase in degrees

phase_rad = phase / 180 *pi