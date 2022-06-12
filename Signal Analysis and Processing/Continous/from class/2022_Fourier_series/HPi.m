function HPi_signal=HPi(t0,T,t)


lower=t0-T/2;
upper=t0+T/2;

if t<lower | t>upper
    HPi_signal=0;
else
    if t==lower | t==upper
        HPi_signal=0.5;
    else
        HPi_signal=1;
    end
end

return
