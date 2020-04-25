function F6 = f6(t,s,ys,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u)

betao_s = beta_is;

if (Ti <= t)
    betao_s = beta_is + ((q1*beta_is -beta_is)/theta1)*(t - Ti);
end
    
if (Ti + theta1 <= t)
    betao_s = q1*beta_is;
end

betai_s = beta_is;

if (Ti <= t)
    betai_s = beta_is + ((u*beta_is -beta_is)/theta)*(t - Ti);
end
    
if (Ti + theta <= t)
    betai_s = u*beta_is;
end

F6 = betai_s*ys*s + betao_s*c_ys*c_s;           % Symptomatic cumulative incidence 