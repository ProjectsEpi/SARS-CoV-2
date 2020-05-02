function G2 = g2(t,e,c_s,c_e,c_ya,c_ys,beta_ia,beta_is,gamma,omega,q1,theta1,Ti)

betao_a = beta_ia;
betao_s = beta_is;

if (Ti <= t)
    betao_a = beta_ia + ((q1*beta_ia -beta_ia)/theta1)*(t - Ti);
    betao_s = beta_is + ((q1*beta_is -beta_is)/theta1)*(t - Ti);
end
    
if (Ti + theta1 <= t)
    betao_a = q1*beta_ia;
    betao_s = q1*beta_is;
end

G2 = (betao_a*c_ya + betao_s*c_ys)*c_s - gamma*c_e + omega*e;