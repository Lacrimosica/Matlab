clc
clear all
close all

N=8;
n=0:N-1;

for h=0:N-1
    for k=0:N-1

        x=exp(1j*k*pi*n/N);
        y=exp(1j*h*pi*n/N);

        z= x * y';
        fprintf("h = %d  k= %d\n" + ...
            "z = %f + %fj\n", h, k, real(z), imag(z));

    end
end

