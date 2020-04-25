function F1 = f11(t,s,ya,ys,N,beta_a,beta_s,rho_t,Lim,Lsm_temp,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
else
    Lsm = Lsm_temp;
end

F1 = -(beta_a*ya + beta_s*ys)*s + rho_t*(Lim/N) + (Lsm/N) - (Lim + Lsm)*(s/N);