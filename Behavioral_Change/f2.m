function F2 = f2(t,s,e,ya,ys,beta_ia,beta_is,gamma,u,omega,theta,Ti)

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

F2 = (beta_a*ya + beta_s*ys)*s - gamma*e - omega*e;