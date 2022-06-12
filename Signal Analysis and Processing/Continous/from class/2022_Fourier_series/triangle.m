function triangle=triangle(t0,T,t)


lower=t0-T;
upper=t0+T;

if t<lower | t>upper
    triangle=0;
else
    if t==lower | t<(lower+T)
        triangle=(t-t0+T)/T;
    else
        triangle=-(t-t0-T)/T;
    end
end

return
