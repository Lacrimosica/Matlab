function exp_smorzato_powt=exp_smorzato_powt(t0,a,pow,f0,t)




if (t-t0)<0
    exp_smorzato_powt=0;
else
    exp_smorzato_powt=exp(-a*(t-t0))*a^(pow+1)*((t-t0).^pow)/gamma(pow+1)*cos(2*pi*f0*(t-t0));
end

return
