% This program solves the mathematical model of the manuscript "Modeling
% behavioral change and COVID-19 containment in Mexico: a trade-off between
% lockdown and compliance" using the 4th order Runge-Kutta method. We employ a csv file 
% which contains a samples of some parameters values (see manuscript for 
% more information) for estimate a confidence band and quantiles of prevalence
% and incidence peaks. The outputs are: 1) prevalence, incidence and cumulative
% incidence and deaths graphs, 2) histograms of the prevalence and incidence peaks
% (values and dates) and 3) files with data of the quantiles (prevalence peaks,
% incidence peaks, cumulative incidence at April 30th and cumulative incidence at August 31st)

clc
clear all

initial_time = datetime('now')

M = csvread('Poisson_samples_v5-5-1.csv',1,0);

S00 = 8918653;  % Total population of Mexico City

%% Parameters for after isolation
omega  = 1/200; % 1/omega is the mean time of application of non-pharmaceutical interventions
q1     = 0.3;   % Desired proportion of reduction of the initial contact rate for out-of-lockdown population
theta  = 150;   % Learning time to reduce contact rate to target level for lockdown population
theta1 = 150;   % Learning time to reduce contact rate to target level for out-of-lockdown population
u      = 0.1;   % Desired proportion of reduction of the initial contact rate for lockdown population
q      = 0.95;  % Fraction of individuals isolated


for id = 1:size(M,1)

    %% Initial conditions
    e0  = 0;
    ya0 = 2;
    ys0 = 2;
    s0  = S00 - ya0 - ys0;
    r0  = 0;
    d0  = 0;
    N0 = S00;

    s  = s0/N0;
    e  = e0/N0;
    ya = ya0/N0;
    ys = ys0/N0;
    r  = r0/N0;

    c_s  = 0;
    c_e  = 0;
    c_ya = 0;
    c_ys = 0;
    c_r  = 0;
    
    d  = d0/N0;
    
    c1 = ya + ys;
    c2 = ys;    

    temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r + d;

    %% Baseline parameters
    beta_ia = M(id,3);     % Contact rate for asymptomatic individuals
    beta_is = M(id,4);     % Contact rate for symptomatic individuals
    gamma   = 1/5.1;       % 1/gamma incubation period
    eta     = M(id,1);     % 1/eta infectious period of asymptomatic infectious individuals
    eta_s   = M(id,2);     % 1/eta infectious period of symptomatic infectious individuals
    rho_i   = 0.45;        % proportion of asymptomatic infectious individuals
    mu      = 0.018;       % mortality rate for disease before Sanitary Emergency Measures

    %% Basic reproductive number
    R0  = ((rho_i*beta_ia/eta) + ((1 - rho_i)*beta_is/eta_s))*s;

    %%
    tim = 0;    % February 17th
    Ti  = 35;   % Date of start isolation - March 23rd
    h   = 0.01; % Step size
    N1  = Ti/h;

    %% Saving initial condition
    XXS(1,id)  = s;             % Susceptible people of the lockdown environment
    XXE(1,id)  = e;             % Exposed people of the lockdown environment
    XXIa(1,id) = ya;            % Asymptomatic people of the lockdown environment
    XXIs(1,id) = ys;            % Symptomatic people of the lockdown environment
    XXR(1,id)  = r;             % Recovered people of the lockdown environment
        
    XXCs(1,id) = c_s;           % Susceptible people of the non-confined environment
    XXCe(1,id) = c_e;           % Exposed people of the non-confined environment
    XXCya(1,id)= c_ya;          % Asymptomatic people of the non-confined environment
    XXCys(1,id)= c_ys;          % Symptomatic people of the non-confined environment
    XXCr(1,id) = c_r;           % Recovered people of the non-confined environment
    
    XXD(1,id)  = d;             % Dead people
    
    XXPrev(1,id) = (ya + ys);   % Prevalence
    
    XXCumu_s(1,id)  = c2;       % Cumulative incidence of symptomatic individuals
    XXCumu_as(1,id) = c1;       % Cumulatice incidence of total (asymptomatic + symptomatic) infected people
    

    %% Solving system for the scenario before the sanitary actions are enforced with 4th order Runge-Kutta method

    for i=1:N1
        
        tim = (i - 1)*h;        

        m1_1 = h*f11(s,ya,ys,beta_ia,beta_is);
        m1_2 = h*f22(s,e,ya,ys,beta_ia,beta_is,gamma);
        m1_3 = h*f33(e,ya,gamma,rho_i,eta);
        m1_4 = h*f44(e,ys,eta_s,rho_i,gamma);
        m1_5 = h*f55(ya,ys,eta,mu,eta_s);
        m1_d = h*f66(ys,mu,eta_s);
        m1_6 = h*g6(tim,s,ya,ys,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m1_7 = h*f6(tim,s,ys,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u);

        m2_1 = h*f11(s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is);
        m2_2 = h*f22(s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,gamma);
        m2_3 = h*f33(e + 0.5*m1_2,ya + 0.5*m1_3,gamma,rho_i,eta);
        m2_4 = h*f44(e + 0.5*m1_2,ys + 0.5*m1_4,eta_s,rho_i,gamma);
        m2_5 = h*f55(ya + 0.5*m1_3,ys + 0.5*m1_4,eta,mu,eta_s);
        m2_d = h*f66(ys + 0.5*m1_4,mu,eta_s);
        m2_6 = h*g6(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m2_7 = h*f6(tim + 0.5*h,s + 0.5*m1_1,ys + 0.5*m1_4,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u);

        m3_1 = h*f11(s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is);
        m3_2 = h*f22(s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,gamma);
        m3_3 = h*f33(e + 0.5*m2_2,ya + 0.5*m2_3,gamma,rho_i,eta);
        m3_4 = h*f44(e + 0.5*m2_2,ys + 0.5*m2_4,eta_s,rho_i,gamma);
        m3_5 = h*f55(ya + 0.5*m2_3,ys + 0.5*m2_4,eta,mu,eta_s);
        m3_d = h*f66(ys + 0.5*m2_4,mu,eta_s);
        m3_6 = h*g6(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m3_7 = h*f6(tim + 0.5*h,s + 0.5*m2_1,ys + 0.5*m2_4,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u);        

        m4_1 = h*f11(s + m3_1,ya + m3_3,ys + m3_4,beta_ia,beta_is);
        m4_2 = h*f22(s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,beta_ia,beta_is,gamma);
        m4_3 = h*f33(e + m3_2,ya + m3_3,gamma,rho_i,eta);
        m4_4 = h*f44(e + m3_2,ys + m3_4,eta_s,rho_i,gamma);
        m4_5 = h*f55(ya + m3_3,ys + m3_4,eta,mu,eta_s);
        m4_d = h*f66(ys + m3_4,mu,eta_s);
        m4_6 = h*g6(tim + h,s + m3_1,ya + m3_3,ys + m3_4,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m4_7 = h*f6(tim + h,s + m3_1,ys + m3_4,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u);          

        s  = s  + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
        e  = e  + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
        ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
        ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
        r  = r  + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);
        d  = d  + (1/6)*(m1_d + 2*m2_d + 2*m3_d + m4_d);
        
        c1 = c1 + (1/6)*(m1_6 + 2*m2_6 + 2*m3_6 + m4_6);
        c2 = c2 + (1/6)*(m1_7 + 2*m2_7 + 2*m3_7 + m4_7);

        temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r + d;

        tim = tim + h;

        %% Saving solution
        if i < N1
                        
            XXS(i+1,id)  = s;             % Susceptible people of the lockdown environment
            XXE(i+1,id)  = e;             % Exposed people of the lockdown environment
            XXIa(i+1,id) = ya;            % Asymptomatic people of the lockdown environment
            XXIs(i+1,id) = ys;            % Symptomatic people of the lockdown environment
            XXR(i+1,id)  = r;             % Recovered people of the lockdown environment

            XXCs(i+1,id) = c_s;           % Susceptible people of the non-confined environment
            XXCe(i+1,id) = c_e;           % Exposed people of the non-confined environment
            XXCya(i+1,id)= c_ya;          % Asymptomatic people of the non-confined environment
            XXCys(i+1,id)= c_ys;          % Symptomatic people of the non-confined environment
            XXCr(i+1,id) = c_r;           % Recovered people of the non-confined environment

            XXD(i+1,id)  = d;             % Dead people

            XXPrev(i+1,id) = (ya + ys);   % Prevalence

            XXCumu_s(i+1,id)  = c2;       % Cumulative incidence of symptomatic individuals
            XXCumu_as(i+1,id) = c1;       % Cumulatice incidence of total (asymptomatic + symptomatic) infected people
       
        end

    end

    %% Initial conditions for dynamics after the sanitary actions are enforced
    temp_s  = q*s;
    temp_e  = q*e;
    temp_ya = q*ya;
    temp_ys = q*ys;
    temp_r  = q*r;

    c_s  = s - temp_s;
    c_e  = e - temp_e;
    c_ya = ya - temp_ya;
    c_ys = ys - temp_ys;
    c_r  = r - temp_r;

    s  = temp_s;
    e  = temp_e;
    ya = temp_ya;
    ys = temp_ys;
    r  = temp_r;
     
    temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r + d;
    
    XXS(N1+1,id)  = s;             % Susceptible people of the lockdown environment
    XXE(N1+1,id)  = e;             % Exposed people of the lockdown environment
    XXIa(N1+1,id) = ya;            % Asymptomatic people of the lockdown environment
    XXIs(N1+1,id) = ys;            % Symptomatic people of the lockdown environment
    XXR(N1+1,id)  = r;             % Recovered people of the lockdown environment
    
    XXCs(N1+1,id) = c_s;           % Susceptible people of the non-confined environment
    XXCe(N1+1,id) = c_e;           % Exposed people of the non-confined environment
    XXCya(N1+1,id)= c_ya;          % Asymptomatic people of the non-confined environment
    XXCys(N1+1,id)= c_ys;          % Symptomatic people of the non-confined environment
    XXCr(N1+1,id) = c_r;           % Recovered people of the non-confined environment

    XXD(N1+1,id)  = d;             % Dead people
    
    XXPrev(N1+1,id) = (ya + ys);   % Prevalence
    
    XXCumu_s(N1+1,id)  = c2;       % Cumulative incidence of symptomatic individuals
    XXCumu_as(N1+1,id) = c1;       % Cumulatice incidence of total (asymptomatic + symptomatic) infected people
    
    j1 = N1+1;
    T  = 196;       % Final time
    N1 = (T - Ti)/h;    
    
    mu = 0.12;                    % mortality rate for disease after Sanitary Emergency Measures
   
    %% Solving system for the scenario after the sanitary actions are enforced with 4th order Runge-Kutta method

    for i=1:N1

        tim = Ti + (i - 1)*h;

        m1_1  = h*f1(tim,s,ya,ys,beta_ia,beta_is,u,omega,theta,Ti);
        m1_2  = h*f2(tim,s,e,ya,ys,beta_ia,beta_is,gamma,u,omega,theta,Ti);
        m1_3  = h*f3(e,ya,gamma,rho_i,eta,omega);
        m1_4  = h*f4(e,ys,eta_s,rho_i,gamma,omega);
        m1_5  = h*f5(ya,ys,r,eta,mu,omega,eta_s);
        m1_6  = h*g1(tim,s,c_s,c_ya,c_ys,beta_ia,beta_is,omega,q1,theta1,Ti);
        m1_7  = h*g2(tim,e,c_s,c_e,c_ya,c_ys,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
        m1_8  = h*g3(ya,c_e,c_ya,gamma,rho_i,eta,omega);
        m1_9  = h*g4(ys,c_e,c_ys,eta_s,rho_i,gamma,omega);
        m1_10 = h*g5(r,c_ya,c_ys,eta,mu,omega,eta_s);
        m1_d  = h*f66(ys + c_ys,mu,eta_s);
        m1_11 = h*g6(tim,s,ya,ys,c_s,c_ya,c_ys,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m1_12 = h*f6(tim,s,ys,c_s,c_ys,beta_is,q1,theta1,Ti,theta,u);

        m2_1  = h*f1(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,u,omega,theta,Ti);
        m2_2  = h*f2(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
        m2_3  = h*f3(e + 0.5*m1_2,ya + 0.5*m1_3,gamma,rho_i,eta,omega);
        m2_4  = h*f4(e + 0.5*m1_2,ys + 0.5*m1_4,eta_s,rho_i,gamma,omega);
        m2_5  = h*f5(ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,eta,mu,omega,eta_s);
        m2_6  = h*g1(tim + 0.5*h,s + 0.5*m1_1,c_s + 0.5*m1_6,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,omega,q1,theta1,Ti);
        m2_7  = h*g2(tim + 0.5*h,e + 0.5*m1_2,c_s + 0.5*m1_6,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
        m2_8  = h*g3(ya + 0.5*m1_3,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,gamma,rho_i,eta,omega);
        m2_9  = h*g4(ys + 0.5*m1_4,c_e + 0.5*m1_7,c_ys + 0.5*m1_9,eta_s,rho_i,gamma,omega);
        m2_10 = h*g5(r + 0.5*m1_5,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,eta,mu,omega,eta_s);
        m2_d  = h*f66(ys + 0.5*m1_4 + c_ys + 0.5*m1_9,mu,eta_s);
        m2_11 = h*g6(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,c_s + 0.5*m1_6,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m2_12 = h*f6(tim + 0.5*h,s + 0.5*m1_1,ys + 0.5*m1_4,c_s + 0.5*m1_6,c_ys + 0.5*m1_9,beta_is,q1,theta1,Ti,theta,u);

        m3_1  = h*f1(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,u,omega,theta,Ti);
        m3_2  = h*f2(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
        m3_3  = h*f3(e + 0.5*m2_2,ya + 0.5*m2_3,gamma,rho_i,eta,omega);
        m3_4  = h*f4(e + 0.5*m2_2,ys + 0.5*m2_4,eta_s,rho_i,gamma,omega);
        m3_5  = h*f5(ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,eta,mu,omega,eta_s);
        m3_6  = h*g1(tim + 0.5*h,s + 0.5*m2_1,c_s + 0.5*m2_6,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,omega,q1,theta1,Ti);
        m3_7  = h*g2(tim + 0.5*h,e + 0.5*m2_2,c_s + 0.5*m2_6,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
        m3_8  = h*g3(ya + 0.5*m2_3,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,gamma,rho_i,eta,omega);
        m3_9  = h*g4(ys + 0.5*m2_4,c_e + 0.5*m2_7,c_ys + 0.5*m2_9,eta_s,rho_i,gamma,omega);
        m3_10 = h*g5(r + 0.5*m2_5,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,eta,mu,omega,eta_s);
        m3_d  = h*f66(ys + 0.5*m2_4 + c_ys + 0.5*m2_9,mu,eta_s);
        m3_11 = h*g6(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,c_s + 0.5*m2_6,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m3_12 = h*f6(tim + 0.5*h,s + 0.5*m2_1,ys + 0.5*m2_4,c_s + 0.5*m2_6,c_ys + 0.5*m2_9,beta_is,q1,theta1,Ti,theta,u);        

        m4_1  = h*f1(tim + h,s + m3_1,ya + m3_3,ys + m3_4,beta_ia,beta_is,u,omega,theta,Ti);
        m4_2  = h*f2(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
        m4_3  = h*f3(e + m3_2,ya + m3_3,gamma,rho_i,eta,omega);
        m4_4  = h*f4(e + m3_2,ys + m3_4,eta_s,rho_i,gamma,omega);
        m4_5  = h*f5(ya + m3_3,ys + m3_4,r + m3_5,eta,mu,omega,eta_s);
        m4_6  = h*g1(tim + h,s + m3_1,c_s + m3_6,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,omega,q1,theta1,Ti);
        m4_7  = h*g2(tim + h,e + m3_2,c_s + m3_6,c_e + m3_7,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
        m4_8  = h*g3(ya + m3_3,c_e + m3_7,c_ya + m3_8,gamma,rho_i,eta,omega);
        m4_9  = h*g4(ys + m3_4,c_e + m3_7,c_ys + m3_9,eta_s,rho_i,gamma,omega);
        m4_10 = h*g5(r + m3_5,c_ya + m3_8,c_ys + m3_9,eta,mu,omega,eta_s);
        m4_d  = h*f66(ys + m3_4 + c_ys + m3_9,mu,eta_s);
        m4_11 = h*g6(tim + h,s + m3_1,ya + m3_3,ys + m3_4,c_s + m3_6,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,q1,theta1,Ti,theta,u);
        m4_12 = h*f6(tim + h,s + m3_1,ys + m3_4,c_s + m3_6,c_ys + m3_9,beta_is,q1,theta1,Ti,theta,u);  


        s  = s  + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
        e  = e  + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
        ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
        ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
        r  = r  + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);

        c_s  = c_s  + (1/6)*(m1_6  + 2*m2_6  + 2*m3_6  + m4_6);
        c_e  = c_e  + (1/6)*(m1_7  + 2*m2_7  + 2*m3_7  + m4_7);
        c_ya = c_ya + (1/6)*(m1_8  + 2*m2_8  + 2*m3_8  + m4_8);
        c_ys = c_ys + (1/6)*(m1_9  + 2*m2_9  + 2*m3_9  + m4_9);
        c_r  = c_r  + (1/6)*(m1_10 + 2*m2_10 + 2*m3_10 + m4_10);
        
        d  = d  + (1/6)*(m1_d + 2*m2_d + 2*m3_d + m4_d);
        
        c1 = c1 + (1/6)*(m1_11 + 2*m2_11 + 2*m3_11 + m4_11);
        c2 = c2 + (1/6)*(m1_12 + 2*m2_12 + 2*m3_12 + m4_12);        

        temp_x = s + e + ya + ys + r + d + c_s + c_e + c_ya + c_ys + c_r;
        
        XXS(j1+1,id)  = s;             % Susceptible people of the lockdown environment
        XXE(j1+1,id)  = e;             % Exposed people of the lockdown environment
        XXIa(j1+1,id) = ya;            % Asymptomatic people of the lockdown environment
        XXIs(j1+1,id) = ys;            % Symptomatic people of the lockdown environment
        XXR(j1+1,id)  = r;             % Recovered people of the lockdown environment
        
        XXCs(j1+1,id) = c_s;           % Susceptible people of the non-confined environment
        XXCe(j1+1,id) = c_e;           % Exposed people of the non-confined environment
        XXCya(j1+1,id)= c_ya;          % Asymptomatic people of the non-confined environment
        XXCys(j1+1,id)= c_ys;          % Symptomatic people of the non-confined environment
        XXCr(j1+1,id) = c_r;           % Recovered people of the non-confined environment
        
        XXD(j1+1,id) = d;              % Dead people of the non-confined environment

        XXPrev(j1+1,id) = (ya + ys + c_ya + c_ys);   % Prevalence

        XXCumu_s(j1+1,id)   = c2;       % Cumulative incidence of symptomatic individuals
        XXCumu_as(j1+1,id)  = c1;       % Cumulatice incidence of total (asymptomatic + symptomatic) infected people
        
        j1 = j1 + 1;        

    end

end

%% Open files for save dates and values of peaks (incidence and prevalence)
file1 = fopen('Cuatiles_dates.csv','w');
file2 = fopen('Cuatiles_value.csv','w');

xx_t = 0:h:T;

j1 = 1;

yy_p(1,:)   = XXPrev(1,:);
yy_is(1,:)  = zeros(1,size(XXCumu_s,2));
yy_ias(1,:) = zeros(1,size(XXCumu_s,2));
yy_cs(1,:)  = XXCumu_s(1,:); 
yy_cas(1,:) = XXCumu_as(1,:);

timo(1)   = datetime(2020,02,17);

%% Data per day
for i=1:size(xx_t,2)
    if xx_t(i) == j1
        yy_p(j1+1,:)    = XXPrev(i,:);
        yy_cs(j1+1,:)   = XXCumu_s(i,:); 
        yy_cas(j1+1,:)  = XXCumu_as(i,:); 
        
        timo(j1+1)    = timo(j1) + caldays(1);
        
        j1       = j1 + 1;
    end
end

%% Calculating incidence values per day
for i = 2:size(yy_cas,1)
    yy_is(i,:)  = yy_cs(i,:) - yy_cs(i-1,:);
    yy_ias(i,:) = yy_cas(i,:) - yy_cas(i-1,:);
end

yy_ia = yy_ias - yy_is;

%% Per day statistics
for i = 1:size(yy_p,1)
    Mat_Prev_days(i,:)    = quantile(yy_p(i,:),[0.025 0.5 0.975]);
    Mat_Inci_s_days(i,:)  = quantile(yy_is(i,:),[0.025 0.5 0.975]);
    Mat_Inci_a_days(i,:)  = quantile(yy_ia(i,:),[0.025 0.5 0.975]);
    Mat_Inci_as_days(i,:) = quantile(yy_ias(i,:),[0.025 0.5 0.975]); 
    Mat_Cumu_s_days(i,:)  = quantile(yy_cs(i,:),[0.025 0.5 0.975]);
    Mat_Cumu_as_days(i,:) = quantile(yy_cas(i,:),[0.025 0.5 0.975]);

end

%% Plotting different characteristics of the disease
figure
subplot(2,3,1); plot(timo, Mat_Prev_days); title('Prevalence')
subplot(2,3,2); plot(timo, Mat_Inci_s_days); title('Incidence - Symptomatic')
subplot(2,3,3); plot(timo, Mat_Inci_as_days); title('Incidence - Both')
subplot(2,3,4); plot(timo, Mat_Cumu_s_days); title('Cumulative - Symptomatic')
subplot(2,3,6); plot(timo, Mat_Cumu_as_days); title('Cumulative - Both')

%% Computing the peak value of incidence and prevalence
for i = 1:size(yy_p,2)
    Peak_prev(i)  = max(yy_p(:,i));        
    Peak_in_s(i)  = max(yy_is(:,i));
    Peak_in_a(i)  = max(yy_ia(:,i));
    Peak_in_as(i) = max(yy_ias(:,i));    
end

quan_peak_prev  = quantile(Peak_prev,[0.025 0.5 0.975]);
quan_peak_in_s  = quantile(Peak_in_s,[0.025 0.5 0.975]);
quan_peak_in_a  = quantile(Peak_in_a,[0.025 0.5 0.975]);
quan_peak_in_as = quantile(Peak_in_as,[0.025 0.5 0.975]);

figure
subplot(1,3,1); hist(Peak_prev); title('Peaks of Prevalence')      % Histogram of the peak values of prevalence
subplot(1,3,2); hist(Peak_in_s); title('Peaks of Symp Incidence')  % Histogram of the peak values of symptomatic incidence
subplot(1,3,3); hist(Peak_in_as); title('Peaks of Incidence')      % Histogram of the peak values of total incidence (asymptomatic + symptomatic)

%% Calculating the dates of peaks (incidence and prevalence)
for i = 1:size(yy_p,2)
    Num_Date_Peak_prev(i)  = find(Peak_prev(i)  == yy_p(:,i));
    Num_Date_Peak_in_s(i)  = find(Peak_in_s(i)  == yy_is(:,i));
    Num_Date_Peak_in_a(i)  = find(Peak_in_a(i)  == yy_ia(:,i));
    Num_Date_Peak_in_as(i) = find(Peak_in_as(i) == yy_ias(:,i));
end

Date_Peak_prev  = timo(Num_Date_Peak_prev(:));
Date_Peak_in_s  = timo(Num_Date_Peak_in_s(:));
Date_Peak_in_a  = timo(Num_Date_Peak_in_a(:));
Date_Peak_in_as = timo(Num_Date_Peak_in_as(:));

quan_peak_prev_date  = round(quantile(Num_Date_Peak_prev,[0.025 0.5 0.975]));      % Calculating the quantiles of dates of the prevalence peaks
quan_peak_in_s_date  = round(quantile(Num_Date_Peak_in_s,[0.025 0.5 0.975]));      % Calculating the quantiles of the dates of the symptomatic incidence peaks
quan_peak_in_a_date  = round(quantile(Num_Date_Peak_in_a,[0.025 0.5 0.975]));      % Calculating the quantiles of the dates of the asymptomatic incidence peaks
quan_peak_in_as_date = round(quantile(Num_Date_Peak_in_as,[0.025 0.5 0.975]));     % Calculating the quantiles of the dates of the total incidence peaks

figure
subplot(1,3,1); histogram(Date_Peak_prev,10);  title('Date of Prevalence Peaks')       % Histogram of the dates of prevalence peak
subplot(1,3,2); histogram(Date_Peak_in_s,10);  title('Date of Symp Incidence Peaks')   % Histogram of the dates of symptomatic incidence peak
subplot(1,3,3); histogram(Date_Peak_in_as,10); title('Date of Incidence Peaks')        % Histogram of the dates of total incidence peak

%% Cumulative incidence at April 30th
cum_is_30  = yy_cs(74,:); 
cum_ias_30 = yy_cas(74,:);

quan_cum_is_30  = quantile(cum_is_30,[0.025 0.5 0.975]); 
quan_cum_ias_30 = quantile(cum_ias_30,[0.025 0.5 0.975]);

quan_cum_is_end  = quantile(yy_cs(end,:),[0.025 0.5 0.975]); 
quan_cum_ias_end = quantile(yy_cas(end,:),[0.025 0.5 0.975]);


%% Saving data of quantiles for prevalence, incidence and cumulative incidence
fprintf(file1,'%s\t %3.0f\t %3.0f\t %3.0f\n','Prevalence',quan_peak_prev_date);
fprintf(file1,'%s\t %3.0f\t %3.0f\t %3.0f\n','Sympto-Inc',quan_peak_in_s_date);
fprintf(file1,'%s\t %3.0f\t %3.0f\t %3.0f\n','Total-Inci',quan_peak_in_as_date);

fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Prevalence',quan_peak_prev);
fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Sympto-Inc',quan_peak_in_s);
fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Total-Inci',quan_peak_in_as);

fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Cumulative-Syntomatic-30',quan_cum_is_30);
fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Cumulative-Total-30',quan_cum_ias_30);
fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Cumulative-Syntomatic-End',quan_cum_is_end);
fprintf(file2,'%s\t %9.9f\t %9.9f\t %9.9f\n','Cumulative-Total-End',quan_cum_ias_end);

period_time = datetime('now') - initial_time
