function rho_t = rho(t,alpha, T0, center, root)
  f0 = (T0)^-1;
  passband = (1-alpha)/(2*f0);
  stopband = (1+alpha)/(2*f0);
  
  rho_t = zeros(1, length(t),'double') + 1i*zeros(1, length(t),'double');
  
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
    endif;
  endif;
endfor;
return

