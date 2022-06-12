function rho_t=rho(t,alpha,T0,center,root)

% t is the array of times at which the signal is evaluated
% alpha is the roll-off factor
% T0 is the 1/2 hieght width
% center is the signal center time
% root at 1 indicates that a root-raised-cosine signal is returned
% root at 0 indicates that a raised-cosine signal is returned

f0=T0^-1;
passband=(1-alpha)/2/f0;
stopband=(1+alpha)/2/f0;

% initializing
rho_t=zeros(1,length(t),'double')+1i*zeros(1,length(t),'double');


for nch=1:length(center)
    
    tt=abs(t-center(nch));
    t=abs(t-center(nch))-passband;

    if alpha==0
        rho_t=(t<=0)+rho_t;
    else
        if root==1
            rho_t=...
            (t<=0)+...
            sqrt(1/2*(1+cos(pi*f0/alpha*(t))).*(t>0).*(abs(tt)<=stopband))+rho_t;
        else
            rho_t=...
            (t<=0)+...
            (1/2*(1+cos(pi*f0/alpha*(t))).*(t>0).*(abs(tt)<=stopband))+rho_t;
        end;
end;

end;

return



