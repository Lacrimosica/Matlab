function gaussiana=gaussiana(t0,a,f0,t)



gaussiana=1/sqrt(2*pi*a^2)*exp(-((t-t0).^2)/2/a^2).*cos(2*pi*f0*(t-t0));

return
