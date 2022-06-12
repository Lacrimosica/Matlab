function triangolo=triangolo_f(t0,T,f0,t)


lower=t0-T;
upper=t0+T;

if t<lower | t>upper
    triangolo=0;
else
    if t==lower | t<(lower+T)
        triangolo=(t-t0+T)/T*cos(2*pi*f0*(t-t0));
    else
        triangolo=-(t-t0-T)/T*cos(2*pi*f0*(t-t0));
    end
end

return
