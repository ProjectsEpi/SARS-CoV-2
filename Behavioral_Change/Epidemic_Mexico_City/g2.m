function G2 = g2(t,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,beta_ia,beta_is,gamma,omega,q1,theta1,Ti,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

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

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

G2 = (beta_a*c_ya + beta_s*c_ys)*c_s - gamma*c_e + omega*e + q_i*(1 - rho_t)*(Lim/N) - (Lim + Lsm)*(c_e/N) - (Lmi + Lms)*Ni*(c_e/N);