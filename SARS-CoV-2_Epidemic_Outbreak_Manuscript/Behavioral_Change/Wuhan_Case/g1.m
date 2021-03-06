function G1 = g1(t,s,c_s,c_ya,c_ys,beta_ia,beta_is,omega,q1,theta1,Ti)

beta_a = beta_ia;
beta_s = beta_is;

if (Ti <= t)
    beta_a = beta_ia + ((q1*beta_ia -beta_ia)/theta1)*(t - Ti);
    beta_s = beta_is + ((q1*beta_is -beta_is)/theta1)*(t - Ti);
end
    
if (Ti + theta1 <= t)
    beta_a = q1*beta_ia;
    beta_s = q1*beta_is;
end

G1 = -(beta_a*c_ya + beta_s*c_ys)*c_s + omega*s;