function F2 = f2(t,s,e,ya,ys,r,beta_ia,beta_is,gamma,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

beta_a = beta_ia;
beta_s = beta_is;

if (Ti <= t)
    beta_a = beta_ia + ((u*beta_ia -beta_ia)/theta)*(t - Ti);
    beta_s = beta_is + ((u*beta_is -beta_is)/theta)*(t - Ti);
end
    
if (Ti + theta <= t)
    beta_a = u*beta_ia;
    beta_s = u*beta_is;
end

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

F2 = (beta_a*ya + beta_s*ys)*s - gamma*e - omega*e - (Lim + Lsm)*(e/N) + (Lmi + Lms)*(1 - Ni)*(e/N);
