function F2 = f22(t,s,e,ya,ys,beta_ia,beta_is,gamma,q1,T_theta,theta)

beta_a = beta_ia;
beta_s = beta_is;

if (T_theta <= t)
    beta_a = beta_ia + ((q1*beta_ia -beta_ia)/theta)*(t - T_theta);
    beta_s = beta_is + ((q1*beta_is -beta_is)/theta)*(t - T_theta);
end
    
if (T_theta + theta <= t)
    beta_a = q1*beta_ia;
    beta_s = q1*beta_is;
end

F2 = (beta_a*ya + beta_s*ys)*s - gamma*e;