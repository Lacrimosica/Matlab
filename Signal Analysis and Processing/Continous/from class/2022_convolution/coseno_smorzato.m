function coseno_smorzato=coseno_smorzato(t0,a,f0,t)




if (t-t0)<0
    coseno_smorzato=0;
else
    coseno_smorzato=exp(-a*(t-t0))*cos(2*pi*f0*(t-t0));
end

return
