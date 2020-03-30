function F2 = f22(t,s,e,ya,ys,N,beta_a,beta_s,gamma,rho_t,q_i,Lim,Lsm_temp,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
else
    Lsm = Lsm_temp;
end

F2 = (beta_a*ya + beta_s*ys)*s + q_i*(1 - rho_t)*(Lim/N) - gamma*e - (Lim + Lsm)*(e/N);