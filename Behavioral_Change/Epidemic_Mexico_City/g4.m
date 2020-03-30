function G4 = g4(t,s,e,ya,ys,r,c_e,c_ys,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

G4 = (1 - rho_i)*gamma*c_e - eta*c_ys + omega*ys - (Lim + Lsm)*(c_ys/N) - (Lmi + Lms)*Ni*(c_ys/N);