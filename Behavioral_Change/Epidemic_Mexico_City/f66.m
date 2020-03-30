function F6 = f66(t,Lim,Lsm_temp,Lmi,Lms_temp,mu,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

F6 = Lim + Lsm - (Lmi + Lms) - mu*N;