clear all
close all
clc

N=10;
n=0:N;
K=3;
P=17;

x = cos(2*pi *n *K /N);
v = cos(2*pi *n *P /N);

figure()
stem(n,x,'o'); hold on
stem(n,v,'s')