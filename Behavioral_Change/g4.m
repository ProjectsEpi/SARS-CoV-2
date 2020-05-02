function G4 = g4(ys,c_e,c_ys,eta,rho_i,gamma,omega)

G4 = (1 - rho_i)*gamma*c_e - eta*c_ys + omega*ys;