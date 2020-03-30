% This program solves system 14-15 of the manuscript "The SARS-CoV-2 epidemic
% outbreak: a review of plausible scenarios of containment and mitigation for
% Mexico" using the 4th order Runge-Kutta method. The outputs are: 
% 1) .txt files with the solution for each instant of time, 2) incidence 
% graph, and 3) cumulative incidence graph for isolated environment.

clc
clear all

initial_time = datetime('now')

%% Initial conditions
s0  = 25210748;
e0  = 0;
ya0 = 0;
ys0 = 0;
r0  = 0;
N0  = s0 + e0 + ya0 + ys0 + r0;

s  = s0/N0;
e  = e0/N0;
ya = ya0/N0;
ys = ys0/N0;
r  = r0/N0;
N  = N0;

temp_x = s + e + ya + ys + r;

%% Baseline parameters
beta_ia = 0.48;         % Contact rate for asymptomatic individuals
beta_is = 0.4;          % Contact rate for symptomatic individuals
gamma   = 1/6;          % 1/gamma incubation period
eta     = 1/10;         % 1/eta infectious period
mu      = 1/(75.1*365); % mortality rate
rho_i   = 0.8;          % proportion of asymptomatic infectious individuals
q_i     = 0.5;          % proportion of air passengers that arrive in the latent stage
N_rho   = 1;            % number of infected international travelers per day

%% Recruitment rates related with the travelling
Lim = 22519;
Lmi = 21838;

Lsm1 = 31276;
Lms1 = 31173;

TN_ms = 34278;
TP_ms = 33168;
TO_ms = 27858;
TS_ms = 14948;

TN_sm = 34318;
TP_sm = 33204;
TO_sm = 27891;
TS_sm = 14975;

Lsm_temp = Lsm1 + TN_sm + TP_sm + TO_sm + TS_sm;
Lms_temp = Lms1 + TN_ms + TP_ms + TO_ms + TS_ms;

h11   = 1.25;               % increase of travelers by holy week
rho_t = 1 - (N_rho/Lim);    % proportion of air infected passengers

%% Basic reproductive number
R01 = sqrt((rho_i*beta_ia + (1 - rho_i)*beta_is)*(s*gamma/((eta + mu)*(gamma + mu))))

%% Open files
file = fopen('dynamics_travelers_non-isolation.txt','w');
file_beta = fopen('transmission_rate.txt','w');

%% 
tim = 0;    % February 21st
Ti  = 31;   % Date of start isolation - March 23rd
h   = 0.01; % Step size
N1  = Ti/h;

%% Saving initial condition and some parameters
fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',tim,s,e,ya,ys,r,temp_x,N);
fprintf(file_beta,'%3.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %9.9f\n',tim,beta_ia,beta_is,beta_ia,beta_is,rho_t,Lim);

%% Solving preisolated System with 4th order Runge-Kutta method

for i=1:N1
    
    tim = (i - 1)*h;
    
    m1_1 = h*f11(tim,s,ya,ys,N,beta_ia,beta_is,rho_t,Lim,Lsm_temp,h11);
    m1_2 = h*f22(tim,s,e,ya,ys,N,beta_ia,beta_is,gamma,rho_t,q_i,Lim,Lsm_temp,h11);
    m1_3 = h*f33(tim,e,ya,N,gamma,rho_i,q_i,rho_t,eta,Lim,Lsm_temp,h11);
    m1_4 = h*f44(tim,e,ys,N,eta,rho_i,gamma,Lim,Lsm_temp,h11);
    m1_5 = h*f55(tim,ya,ys,r,N,eta,Lim,Lsm_temp,h11);
    m1_6 = h*f66(tim,Lim,Lsm_temp,Lmi,Lms_temp,mu,N,h11);
    
    m2_1 = h*f11(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,N + 0.5*m1_6,beta_ia,beta_is,rho_t,Lim,Lsm_temp,h11);
    m2_2 = h*f22(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,N + 0.5*m1_6,beta_ia,beta_is,gamma,rho_t,q_i,Lim,Lsm_temp,h11);
    m2_3 = h*f33(tim + 0.5*h,e + 0.5*m1_2,ya + 0.5*m1_3,N + 0.5*m1_6,gamma,rho_i,q_i,rho_t,eta,Lim,Lsm_temp,h11);
    m2_4 = h*f44(tim + 0.5*h,e + 0.5*m1_2,ys + 0.5*m1_4,N + 0.5*m1_6,eta,rho_i,gamma,Lim,Lsm_temp,h11);
    m2_5 = h*f55(tim + 0.5*h,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,N + 0.5*m1_6,eta,Lim,Lsm_temp,h11);
    m2_6 = h*f66(tim + 0.5*h,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + 0.5*m1_6,h11);    
    
    m3_1 = h*f11(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,N + 0.5*m2_6,beta_ia,beta_is,rho_t,Lim,Lsm_temp,h11);
    m3_2 = h*f22(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,N + 0.5*m2_6,beta_ia,beta_is,gamma,rho_t,q_i,Lim,Lsm_temp,h11);
    m3_3 = h*f33(tim + 0.5*h,e + 0.5*m2_2,ya + 0.5*m2_3,N + 0.5*m2_6,gamma,rho_i,q_i,rho_t,eta,Lim,Lsm_temp,h11);
    m3_4 = h*f44(tim + 0.5*h,e + 0.5*m2_2,ys + 0.5*m2_4,N + 0.5*m2_6,eta,rho_i,gamma,Lim,Lsm_temp,h11);
    m3_5 = h*f55(tim + 0.5*h,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,N + 0.5*m2_6,eta,Lim,Lsm_temp,h11);
    m3_6 = h*f66(tim + 0.5*h,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + 0.5*m2_6,h11);    
    
    m4_1 = h*f11(tim + h,s + m3_1,ya + m3_3,ys + m3_4,N + m3_6,beta_ia,beta_is,rho_t,Lim,Lsm_temp,h11);
    m4_2 = h*f22(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,N + m3_6,beta_ia,beta_is,gamma,rho_t,q_i,Lim,Lsm_temp,h11);
    m4_3 = h*f33(tim + h,e + m3_2,ya + m3_3,N + m3_6,gamma,rho_i,q_i,rho_t,eta,Lim,Lsm_temp,h11);
    m4_4 = h*f44(tim + h,e + m3_2,ys + m3_4,N + m3_6,eta,rho_i,gamma,Lim,Lsm_temp,h11);
    m4_5 = h*f55(tim + h,ya + m3_3,ys + m3_4,r + m3_5,N + m3_6,eta,Lim,Lsm_temp,h11);
    m4_6 = h*f66(tim + h,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + m3_6,h11);   
           
        
    s  = s  + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
    e  = e  + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
    ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
    ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
    r  = r  + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);
    N  = N  + (1/6)*(m1_6 + 2*m2_6 + 2*m3_6 + m4_6);
    
    temp_x = s + e + ya + ys + r;
    
 
    %% Saving solution and some parameters
    fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',tim,s,e,ya,ys,r,temp_x,N);
    fprintf(file_beta,'%3.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %9.9f\n',tim,beta_ia,beta_is,beta_ia,beta_is,rho_t,Lim);

end

fclose(file);

%% Parameters for after isolation
omega  = 1/100;     % 1/omega is the mean time of application of non-pharmaceutical interventions
q1     = 0.5;       % Desired proportion of reduction of the initial contact rate for non-isolated environment
theta  = 7;        % Learning time to reduce contact rate to target level for isolated environment
theta1 = 15;        % Learning time to reduce contact rate to target level for non-isolated environment
u      = 0.4;      % Desired proportion of reduction of the initial contact rate for isolated environment
q      = 0.8;       % Fraction of individuals in the population goes into isolation.

T  = 365;
N1 = (T-Ti)/h;

%% Travelers parameters for dynamics after isolation
h22 = 0.01;     % decrease of travelers by isolation
h33 = 0.5;      % reduction of infected international travelers

Lim = h22*22519;
Lmi = h22*21838;

Lsm1 = 31276;
Lms1 = 31173;

TN_ms = 34278;
TP_ms = 33168;
TO_ms = 27858;
TS_ms = 14948;

T_temp_ms = TN_ms + TP_ms + TO_ms + TS_ms;

TN_sm = 34318;
TP_sm = 33204;
TO_sm = 27891;
TS_sm = 14975;

Lsm_temp = h22*(Lsm1 + TN_sm + TP_sm + TO_sm + TS_sm);
Lms_temp = h22*(Lms1 + TN_ms + TP_ms + TO_ms + TS_ms);

h11 = 1.25;                     % increase of travelers by holy week
rho_t = 1 - (h33*N_rho/Lim);    % proportion of air infected passengers after isolation

%% Initial conditions for dynamics with isolation
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

temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r;

%% Open files
file_with = fopen('dynamics_travelers_isolation.txt','w');

%% Saving data
fprintf(file_with,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',tim,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,N);

%% Solving dynamics with isolation

for i=1:N1

    tim = Ti + (i - 1)*h;

    m1_1  = h*f1(tim,s,e,ya,ys,r,beta_ia,beta_is,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_2  = h*f2(tim,s,e,ya,ys,r,beta_ia,beta_is,gamma,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_3  = h*f3(tim,s,e,ya,ys,r,gamma,rho_i,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_4  = h*f4(tim,s,e,ya,ys,r,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_5  = h*f5(tim,s,e,ya,ys,r,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_6  = h*g1(tim,s,e,ya,ys,r,c_s,c_ya,c_ys,beta_ia,beta_is,omega,q1,theta1,Ti,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_7  = h*g2(tim,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,beta_ia,beta_is,gamma,omega,q1,theta1,Ti,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_8  = h*g3(tim,s,e,ya,ys,r,c_e,c_ya,gamma,rho_i,eta,omega,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_9  = h*g4(tim,s,e,ya,ys,r,c_e,c_ys,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_10 = h*g5(tim,s,e,ya,ys,r,c_ya,c_ys,c_r,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11);
    m1_N  = h*f6(tim,s,e,ya,ys,r,Lim,Lsm_temp,Lmi,Lms_temp,mu,N,h11);

    m2_1  = h*f1(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,beta_ia,beta_is,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_2  = h*f2(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,beta_ia,beta_is,gamma,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_3  = h*f3(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,gamma,rho_i,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_4  = h*f4(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_5  = h*f5(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_6  = h*g1(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,c_s + 0.5*m1_6,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,omega,q1,theta1,Ti,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_7  = h*g2(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,c_s + 0.5*m1_6,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_8  = h*g3(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,gamma,rho_i,eta,omega,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_9  = h*g4(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,c_e + 0.5*m1_7,c_ys + 0.5*m1_9,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_10 = h*g5(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,c_r + 0.5*m1_10,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m1_N,h11);
    m2_N  = h*f6(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + 0.5*m1_N,h11);

    m3_1  = h*f1(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,beta_ia,beta_is,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_2  = h*f2(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,beta_ia,beta_is,gamma,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_3  = h*f3(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,gamma,rho_i,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_4  = h*f4(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_5  = h*f5(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_6  = h*g1(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,c_s + 0.5*m2_6,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,omega,q1,theta1,Ti,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_7  = h*g2(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,c_s + 0.5*m2_6,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_8  = h*g3(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,gamma,rho_i,eta,omega,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_9  = h*g4(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,c_e + 0.5*m2_7,c_ys + 0.5*m2_9,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_10 = h*g5(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,c_r + 0.5*m2_10,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N+ 0.5*m2_N,h11);
    m3_N  = h*f6(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + 0.5*m2_N,h11);

    m4_1  = h*f1(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,beta_ia,beta_is,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_2  = h*f2(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,beta_ia,beta_is,gamma,u,omega,theta,Ti,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_3  = h*f3(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,gamma,rho_i,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_4  = h*f4(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_5  = h*f5(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_6  = h*g1(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,c_s + m3_6,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,omega,q1,theta1,Ti,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_7  = h*g2(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,c_s + m3_6,c_e + m3_7,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_8  = h*g3(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,c_e + m3_7,c_ya + m3_8,gamma,rho_i,eta,omega,q_i,rho_t,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_9  = h*g4(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,c_e + m3_7,c_ys + m3_9,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_10 = h*g5(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,c_ya + m3_8,c_ys + m3_9,c_r + m3_10,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N + m3_N,h11);
    m4_N  = h*f6(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,r + m3_5,Lim,Lsm_temp,Lmi,Lms_temp,mu,N + m3_N,h11);

    s  = s + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
    e  = e + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
    ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
    ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
    r  = r + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);

    c_s  = c_s + (1/6)*(m1_6 + 2*m2_6 + 2*m3_6 + m4_6);
    c_e  = c_e + (1/6)*(m1_7 + 2*m2_7 + 2*m3_7 + m4_7);
    c_ya = c_ya + (1/6)*(m1_8 + 2*m2_8 + 2*m3_8 + m4_8);
    c_ys = c_ys + (1/6)*(m1_9 + 2*m2_9 + 2*m3_9 + m4_9);
    c_r  = c_r + (1/6)*(m1_10 + 2*m2_10 + 2*m3_10 + m4_10);
    
    N = N + (1/6)*(m1_N + 2*m2_N + 2*m3_N + m4_N);

    temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r;
    
    %% Saving solutions
    fprintf(file_with,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',(tim + h),s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,N);
            
    %% Time-dependent contact rates
    if (Ti <= (Ti + i*h))
        beta_a1 = beta_ia + ((u*beta_ia -beta_ia)/theta)*((Ti+i*h) - Ti);
        beta_s1 = beta_is + ((u*beta_is -beta_is)/theta)*((Ti+i*h) - Ti);
                
        beta_a2 = beta_ia + ((q1*beta_ia -beta_ia)/theta1)*((Ti+i*h) - Ti);
        beta_s2 = beta_is + ((q1*beta_is -beta_is)/theta1)*((Ti+i*h) - Ti);                
    end

    if (Ti + theta <= (Ti + i*h))       
        beta_a1 = u*beta_ia;
        beta_s1 = u*beta_is;
    end

    if (Ti + theta1 <= (Ti + i*h))
        beta_a2 = q1*beta_ia;
        beta_s2 = q1*beta_is;
    end
    
    %% Saving some parameters
    fprintf(file_beta,'%3.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\t %9.9f\n',(tim + h),beta_a1,beta_s1,beta_a2,beta_s2,rho_t,Lim);

end
        
fclose(file_with);
fclose(file_beta);

%% Incidence Graphic

delimiterIn = ' ';

file_non  = 'dynamics_travelers_non-isolation.txt';
file_with = 'dynamics_travelers_isolation.txt';
file_beta =  'transmission_rate.txt';

A_noni = importdata(file_non,delimiterIn);
A_with = importdata(file_with,delimiterIn);
A_beta = importdata(file_beta,delimiterIn);

figure
plot(A_noni(1:100*Ti + 1,1),A_noni(1:100*Ti + 1,end).*((A_beta(1:100*Ti + 1,2).*A_noni(1:100*Ti + 1,4) + A_beta(1:100*Ti + 1,3).*A_noni(1:100*Ti + 1,5)).*A_noni(1:100*Ti + 1,2) + (1 - A_beta(1:100*Ti + 1,end - 1)).*(A_beta(1:100*Ti + 1,end)./A_noni(1:100*Ti + 1,end))),'b')
hold on
plot(A_with(:,1),A_with(:,end).*(A_beta(100*Ti + 1:end,2).*A_with(:,4) + A_beta(100*Ti + 1:end,3).*A_with(:,5)).*A_with(:,2),'r--')   % Total infectious people (proportion)
ylabel('Incidence')
xlabel('Time')

%% Cumulative incidence for isolated enviroment
j1 = 1;
for i=1:size(A_noni,1)
    if A_noni(i,1) == j1
        temp_non(j1) = A_noni(i,end)*((A_beta(i,2).*A_noni(i,4) + A_beta(i,3)*A_noni(i,5))*A_noni(i,2) + (1 - A_beta(i,end - 1))*(A_beta(i,end)/A_noni(i,end)));
        j1 = j1 + 1;
    end
end

j1 = Ti + 1;
j2 = 1;
for i=1:size(A_with,1)
    if A_with(i,1) == j1
        temp_with_1(j2) = A_with(i,end).*(A_beta(i,2).*A_with(i,4) + A_beta(i,3).*A_with(i,5)).*A_with(i,2);
        j1 = j1 + 1;
        j2 = j2 + 1;
    end
end

final_cum_1 = [temp_non temp_with_1]; 
final_time = 1:size(final_cum_1,2);

temp_v1 = 0;

for i=1:size(final_time,2)
    temp_v1 = temp_v1 + final_cum_1(i);
    final_cum_in(i)  = temp_v1;
end

figure
plot(final_time,final_cum_in,'r')
title('Cumulative Incidence')
xlabel('Time')

period_time = datetime('now') - initial_time