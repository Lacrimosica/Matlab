function porta=porta_f(t0,T,f0,t)


lower=t0-T/2;
upper=t0+T/2;

if t<lower | t>upper
    porta=0;
else
    if t==lower | t==upper
        porta=0.5*cos(2*pi*f0*(t-t0));
    else
        porta=1*cos(2*pi*f0*(t-t0));
    end
end

return
