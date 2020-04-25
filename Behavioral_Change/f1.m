function F1 = f1(t,s,ya,ys,beta_ia,beta_is,u,omega,theta,Ti)

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

F1 = -(beta_a*ya + beta_s*ys)*s - omega*s;