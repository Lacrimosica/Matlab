function exp_smorzato=exp_smorzato(t0,a,t)




if (t-t0)<0
    exp_smorzato=0;
else
    exp_smorzato=exp(-a*(t-t0));
end

return
