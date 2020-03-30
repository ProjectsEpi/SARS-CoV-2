% This program solves system 8 of the manuscript "The SARS-CoV-2 epidemic
% outbreak: a review of plausible scenarios of containment and mitigation for
% Mexico" using the 4th order Runge-Kutta method. The outputs are: 
% 1) .txt files with the solution for each instant of time, 2) prevalence graph,
% and 3) cumulative incidence graph.

clc
clear all

initial_time = datetime('now')

%% Initial conditions
s0  = 25210738;
e0  = 0;
ya0 = 10;
ys0 = 0;
r0  = 0;
N0 = s0 + ya0;

s  = s0/N0;
e  = e0/N0;
ya = ya0/N0;
ys = ys0/N0;
r  = r0/N0;

temp_x = s + e + ya + ys + r;

%% Baseline parameters
beta_ia = 0.48;         % Contact rate for asymptomatic individuals
beta_is = 0.4;          % Contact rate for symptomatic individuals
gamma   = 1/6;          % 1/gamma incubation period
eta     = 1/10;         % 1/eta infectious period
rho_i   = 0.8;          % proportion of asymptomatic infectious individuals

%% Parameters for contact rates
q1      = 0.5;          % Desired proportion of reduction of the initial contact rate
T_theta = 15;           % Time delay for the application of contact rate reduction
theta   = 180;          % Learning time to reduce contact rate to target level

%% Open files
file = fopen('dynamics_behavioral_change.txt','w');

%%
T    = 250;
tim  = 0;    % February 21st
h    = 0.01; % Step size
N1   = T/h;

%% Saving data
fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\n',tim,s,e,ya,ys,r,beta_ia,beta_is);

%% Solving EDO System with 4th order Runge-Kutta method

for i=1:N1
    
    tim = (i - 1)*h;
    
    m1_1 = h*f11(tim,s,ya,ys,beta_ia,beta_is,q1,T_theta,theta);
    m1_2 = h*f22(tim,s,e,ya,ys,beta_ia,beta_is,gamma,q1,T_theta,theta);
    m1_3 = h*f33(e,ya,gamma,rho_i,eta);
    m1_4 = h*f44(e,ys,eta,rho_i,gamma);
    m1_5 = h*f55(ya,ys,eta);

    m2_1 = h*f11(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,q1,T_theta,theta);
    m2_2 = h*f22(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,gamma,q1,T_theta,theta);
    m2_3 = h*f33(e + 0.5*m1_2,ya + 0.5*m1_3,gamma,rho_i,eta);
    m2_4 = h*f44(e + 0.5*m1_2,ys + 0.5*m1_4,eta,rho_i,gamma);
    m2_5 = h*f55(ya + 0.5*m1_3,ys + 0.5*m1_4,eta);
    
    m3_1 = h*f11(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,q1,T_theta,theta);
    m3_2 = h*f22(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,gamma,q1,T_theta,theta);
    m3_3 = h*f33(e + 0.5*m2_2,ya + 0.5*m2_3,gamma,rho_i,eta);
    m3_4 = h*f44(e + 0.5*m2_2,ys + 0.5*m2_4,eta,rho_i,gamma);
    m3_5 = h*f55(ya + 0.5*m2_3,ys + 0.5*m2_4,eta);   
    
    m4_1 = h*f11(tim + h,s + m3_1,ya + m3_3,ys + m3_4,beta_ia,beta_is,q1,T_theta,theta);
    m4_2 = h*f22(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,beta_ia,beta_is,gamma,q1,T_theta,theta);
    m4_3 = h*f33(e + m3_2,ya + m3_3,gamma,rho_i,eta);
    m4_4 = h*f44(e + m3_2,ys + m3_4,eta,rho_i,gamma);
    m4_5 = h*f55(ya + m3_3,ys + m3_4,eta);  
    
    s  = s  + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
    e  = e  + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
    ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
    ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
    r  = r  + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);

    temp_x = s + e + ya + ys + r;
        
    %% Time-dependent contact rates
    
    beta_fa = beta_ia;
    beta_fs = beta_is;    
    
    if (T_theta <= (tim + h))
        beta_Ta = q1*beta_ia;
        beta_fa = beta_ia + ((beta_Ta -beta_ia)/theta)*(tim + h - T_theta);
        
        beta_Ts = q1*beta_is;
        beta_fs = beta_is + ((beta_Ts -beta_is)/theta)*(tim + h - T_theta);
    end
    
    if (T_theta + theta <= (tim + h))
        beta_fa = q1*beta_ia;
        beta_fs = q1*beta_is;
    end
    
    
    %% Saving solution
    fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\n',tim + h,s,e,ya,ys,r,beta_fa,beta_fs);    

end

%% Close files
fclose(file);

%% Prevalence Graphic
delimiterIn = ' ';

file = 'dynamics_behavioral_change.txt';
A    = importdata(file,delimiterIn);

figure
plot(A(:,1),A(:,4) + A(:,5),'b')

j1 = 1;
j_temp = 1;

%% Cumulative incidence
for i=1:size(A,1)
    if A(i,1) == j1
        temp     = (A(i,end-1)*A(i,4) + A(i,end)*A(i,5))*A(i,2);
        sum2(j1) = temp;
        j_temp   = i + 1;
        j1       = j1 + 1;
    end
end

sum1(1) = sum2(1);
v(1)    = 1;

for i=2:size(sum2,2)
    sum1(i) = sum1(i - 1) + sum2(i);
    v(i)    = i;
end

figure
plot(v,sum1,'r')

period_time = datetime('now') - initial_time