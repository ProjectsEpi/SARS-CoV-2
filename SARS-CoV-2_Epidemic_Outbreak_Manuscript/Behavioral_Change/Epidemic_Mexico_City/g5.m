function G5 = g5(t,s,e,ya,ys,r,c_ya,c_ys,c_r,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

G5 = eta*(c_ya + c_ys) + omega*r - (Lim + Lsm)*(c_r/N) - (Lmi + Lms)*Ni*(c_r/N);