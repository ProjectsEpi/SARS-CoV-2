function F5 = f55(t,ya,ys,r,N,eta,Lim,Lsm_temp,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
else
    Lsm = Lsm_temp;
end

F5 = eta*(ya + ys) - (Lim + Lsm)*(r/N);