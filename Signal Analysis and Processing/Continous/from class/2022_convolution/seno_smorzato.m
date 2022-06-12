function seno_smorzato=seno_smorzato(t0,a,f0,t)




if (t-t0)<0
    seno_smorzato=0;
else
    seno_smorzato=exp(-a*(t-t0))*sin(2*pi*f0*(t-t0));
end

return
