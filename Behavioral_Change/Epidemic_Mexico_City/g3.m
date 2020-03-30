function G3 = g3(t,s,e,ya,ys,r,c_e,c_ya,gamma,rho_i,eta,omega,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

G3 = rho_i*gamma*c_e - eta*c_ya + omega*ya + (1 - q_i)*(1 - rho_t)*(Lim/N) - (Lim + Lsm)*(c_ya/N) - (Lmi + Lms)*Ni*(c_ya/N);