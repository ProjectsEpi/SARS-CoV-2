function G6 = g6(t,s,ya,ys,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u)

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

betai_a = beta_ia;
betai_s = beta_is;

if (Ti <= t)
    betai_a = beta_ia + ((u*beta_ia -beta_ia)/theta)*(t - Ti);
    betai_s = beta_is + ((u*beta_is -beta_is)/theta)*(t - Ti);
end
    
if (Ti + theta <= t)
    betai_a = u*beta_ia;
    betai_s = u*beta_is;
end

G6 = (betai_a*ya + betai_s*ys)*s + (betao_a*c_ya + betao_s*c_ys)*c_s;           % Total cumulative incidence 