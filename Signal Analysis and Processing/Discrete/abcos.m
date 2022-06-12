clear all
close all

clc
N=10;
A=2;
B=3;

n=0:2*N-1;

x = A*cos(2*pi*n/N) + B*sin(4*pi*n/N);

stem(n,x); grid on