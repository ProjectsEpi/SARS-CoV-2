% This program solves system 10-12 of the manuscript "The SARS-CoV-2 epidemic
% outbreak: a review of plausible scenarios of containment and mitigation for
% Mexico" using the 4th order Runge-Kutta method. The outputs are: 
% 1) .txt files with the solution for each instant of time, 2) prevalence graph,
% 3) incidence graph, and 4) cumulative incidence vs Hubei JHU data graph
% for isolated environment.


clc
clear all

initial_time = datetime('now')

%% Initial conditions
s0  = 0.999994;
e0  = 0;
ya0 = 1 - s0;
ys0 = 0;
r0  = 0;

c_s0  = 0;
c_e0  = 0;
c_ya0 = 0;
c_ys0 = 0;
c_r0  = 0;

s  = s0;
e  = e0;
ya = ya0;
ys = ys0;
r  = r0;

c_s  = c_s0;
c_e  = c_e0;
c_ya = c_ya0;
c_ys = c_ys0;
c_r  = c_r0;

temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r;

%% Baseline parameters
beta_ia = 0.65;     % Initial contact rate for asymptomatic individuals
beta_is = 0.55;     % Initial contact rate for symptomatic individuals
gamma   = 1/6;      % 1/gamma incubation period
eta     = 1/10;     % 1/eta infectious period
rho_i   = 0.8;      % proportion of asymptomatic infectious individuals

%% Basic reproductive number
R0  = sqrt((rho_i*beta_ia + (1 - rho_i)*beta_is)*(s0/(eta)))

%% Open files
file = fopen('dynamics_nonisolation.txt','w');

%%
tim = 0;    % December 31st
Ti  = 26;   % Date of start isolation - January 26th
h   = 0.01; % Step size
N1  = Ti/h;

%% Saving initial condition
fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.6f\t %1.6f\n',tim,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,beta_ia,beta_is);

%% Solving preisolated System with 4th order Runge-Kutta method

for i=1:N1
    
    m1_1 = h*f11(s,ya,ys,beta_ia,beta_is);
    m1_2 = h*f22(s,e,ya,ys,beta_ia,beta_is,gamma);
    m1_3 = h*f33(e,ya,gamma,rho_i,eta);
    m1_4 = h*f44(e,ys,eta,rho_i,gamma);
    m1_5 = h*f55(ya,ys,eta);
    
    m2_1 = h*f11(s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is);
    m2_2 = h*f22(s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,gamma);
    m2_3 = h*f33(e + 0.5*m1_2,ya + 0.5*m1_3,gamma,rho_i,eta);
    m2_4 = h*f44(e + 0.5*m1_2,ys + 0.5*m1_4,eta,rho_i,gamma);
    m2_5 = h*f55(ya + 0.5*m1_3,ys + 0.5*m1_4,eta);    
        
    m3_1 = h*f11(s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is);
    m3_2 = h*f22(s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,gamma);
    m3_3 = h*f33(e + 0.5*m2_2,ya + 0.5*m2_3,gamma,rho_i,eta);
    m3_4 = h*f44(e + 0.5*m2_2,ys + 0.5*m2_4,eta,rho_i,gamma);
    m3_5 = h*f55(ya + 0.5*m2_3,ys + 0.5*m2_4,eta); 
        
    m4_1 = h*f11(s + m3_1,ya + m3_3,ys + m3_4,beta_ia,beta_is);
    m4_2 = h*f22(s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,beta_ia,beta_is,gamma);
    m4_3 = h*f33(e + m3_2,ya + m3_3,gamma,rho_i,eta);
    m4_4 = h*f44(e + m3_2,ys + m3_4,eta,rho_i,gamma);
    m4_5 = h*f55(ya + m3_3,ys + m3_4,eta); 
        
    s  = s  + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
    e  = e  + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
    ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
    ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
    r  = r  + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);

    temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r;
    
    tim = tim + h;
    
    %% Saving solution
    fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.6f\t %1.6f\n',tim,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,beta_ia,beta_is);
       
end

%% Close files
fclose(file);

%% Parameters for after isolation
omega  = 0.01;  % 1/omega is the mean time of application of non-pharmaceutical interventions
q1     = 0.5;   % Desired proportion of reduction of the initial contact rate for non-isolated environment
theta  = 7;     % Learning time to reduce contact rate to target level for isolated environment
theta1 = 15;    % Learning time to reduce contact rate to target level for non-isolated environment
u      = 0.1;   % Desired proportion of reduction of the initial contact rate for isolated environment
q      = 0.99;  % Fraction of individuals in the population goes into isolation.

T  = 120;       % Final time
N1 = (T - Ti)/h;


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
        
file_with = fopen('dynamics_withisolation.txt','w');
fprintf(file_with,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\n',tim,s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,beta_ia,beta_is,beta_ia,beta_is);

%% Solving dynamics with isolation

for i=1:N1
    
    tim = Ti + (i - 1)*h;

    m1_1  = h*f1(tim,s,ya,ys,beta_ia,beta_is,u,omega,theta,Ti);
    m1_2  = h*f2(tim,s,e,ya,ys,beta_ia,beta_is,gamma,u,omega,theta,Ti);
    m1_3  = h*f3(e,ya,gamma,rho_i,eta,omega);
    m1_4  = h*f4(e,ys,eta,rho_i,gamma,omega);
    m1_5  = h*f5(ya,ys,r,eta,omega);
    m1_6  = h*g1(tim,s,c_s,c_ya,c_ys,beta_ia,beta_is,omega,q1,theta1,Ti);
    m1_7  = h*g2(tim,e,c_s,c_e,c_ya,c_ys,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
    m1_8  = h*g3(ya,c_e,c_ya,gamma,rho_i,eta,omega);
    m1_9  = h*g4(ys,c_e,c_ys,eta,rho_i,gamma,omega);
    m1_10 = h*g5(r,c_ya,c_ys,eta,omega);

    m2_1  = h*f1(tim + 0.5*h,s + 0.5*m1_1,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,u,omega,theta,Ti);
    m2_2  = h*f2(tim + 0.5*h,s + 0.5*m1_1,e + 0.5*m1_2,ya + 0.5*m1_3,ys + 0.5*m1_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
    m2_3  = h*f3(e + 0.5*m1_2,ya + 0.5*m1_3,gamma,rho_i,eta,omega);
    m2_4  = h*f4(e + 0.5*m1_2,ys + 0.5*m1_4,eta,rho_i,gamma,omega);
    m2_5  = h*f5(ya + 0.5*m1_3,ys + 0.5*m1_4,r + 0.5*m1_5,eta,omega);
    m2_6  = h*g1(tim + 0.5*h,s + 0.5*m1_1,c_s + 0.5*m1_6,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,omega,q1,theta1,Ti);
    m2_7  = h*g2(tim + 0.5*h,e + 0.5*m1_2,c_s + 0.5*m1_6,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
    m2_8  = h*g3(ya + 0.5*m1_3,c_e + 0.5*m1_7,c_ya + 0.5*m1_8,gamma,rho_i,eta,omega);
    m2_9  = h*g4(ys + 0.5*m1_4,c_e + 0.5*m1_7,c_ys + 0.5*m1_9,eta,rho_i,gamma,omega);
    m2_10 = h*g5(r + 0.5*m1_5,c_ya + 0.5*m1_8,c_ys + 0.5*m1_9,eta,omega);

    m3_1  = h*f1(tim + 0.5*h,s + 0.5*m2_1,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,u,omega,theta,Ti);
    m3_2  = h*f2(tim + 0.5*h,s + 0.5*m2_1,e + 0.5*m2_2,ya + 0.5*m2_3,ys + 0.5*m2_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
    m3_3  = h*f3(e + 0.5*m2_2,ya + 0.5*m2_3,gamma,rho_i,eta,omega);
    m3_4  = h*f4(e + 0.5*m2_2,ys + 0.5*m2_4,eta,rho_i,gamma,omega);
    m3_5  = h*f5(ya + 0.5*m2_3,ys + 0.5*m2_4,r + 0.5*m2_5,eta,omega);
    m3_6  = h*g1(tim + 0.5*h,s + 0.5*m2_1,c_s + 0.5*m2_6,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,omega,q1,theta1,Ti);
    m3_7  = h*g2(tim + 0.5*h,e + 0.5*m2_2,c_s + 0.5*m2_6,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
    m3_8  = h*g3(ya + 0.5*m2_3,c_e + 0.5*m2_7,c_ya + 0.5*m2_8,gamma,rho_i,eta,omega);
    m3_9  = h*g4(ys + 0.5*m2_4,c_e + 0.5*m2_7,c_ys + 0.5*m2_9,eta,rho_i,gamma,omega);
    m3_10 = h*g5(r + 0.5*m2_5,c_ya + 0.5*m2_8,c_ys + 0.5*m2_9,eta,omega);

    m4_1  = h*f1(tim + h,s + m3_1,ya + m3_3,ys + m3_4,beta_ia,beta_is,u,omega,theta,Ti);
    m4_2  = h*f2(tim + h,s + m3_1,e + m3_2,ya + m3_3,ys + m3_4,beta_ia,beta_is,gamma,u,omega,theta,Ti);
    m4_3  = h*f3(e + m3_2,ya + m3_3,gamma,rho_i,eta,omega);
    m4_4  = h*f4(e + m3_2,ys + m3_4,eta,rho_i,gamma,omega);
    m4_5  = h*f5(ya + m3_3,ys + m3_4,r + m3_5,eta,omega);
    m4_6  = h*g1(tim + h,s + m3_1,c_s + m3_6,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,omega,q1,theta1,Ti);
    m4_7  = h*g2(tim + h,e + m3_2,c_s + m3_6,c_e + m3_7,c_ya + m3_8,c_ys + m3_9,beta_ia,beta_is,gamma,omega,q1,theta1,Ti);
    m4_8  = h*g3(ya + m3_3,c_e + m3_7,c_ya + m3_8,gamma,rho_i,eta,omega);
    m4_9  = h*g4(ys + m3_4,c_e + m3_7,c_ys + m3_9,eta,rho_i,gamma,omega);
    m4_10 = h*g5(r + m3_5,c_ya + m3_8,c_ys + m3_9,eta,omega);


    s  = s + (1/6)*(m1_1 + 2*m2_1 + 2*m3_1 + m4_1);
    e  = e + (1/6)*(m1_2 + 2*m2_2 + 2*m3_2 + m4_2);
    ya = ya + (1/6)*(m1_3 + 2*m2_3 + 2*m3_3 + m4_3);
    ys = ys + (1/6)*(m1_4 + 2*m2_4 + 2*m3_4 + m4_4);
    r  = r + (1/6)*(m1_5 + 2*m2_5 + 2*m3_5 + m4_5);

    c_s = c_s + (1/6)*(m1_6 + 2*m2_6 + 2*m3_6 + m4_6);
    c_e = c_e + (1/6)*(m1_7 + 2*m2_7 + 2*m3_7 + m4_7);
    c_ya = c_ya + (1/6)*(m1_8 + 2*m2_8 + 2*m3_8 + m4_8);
    c_ys = c_ys + (1/6)*(m1_9 + 2*m2_9 + 2*m3_9 + m4_9);
    c_r = c_r + (1/6)*(m1_10 + 2*m2_10 + 2*m3_10 + m4_10);

    temp_x = s + e + ya + ys + r + c_s + c_e + c_ya + c_ys + c_r;

    %% Time-dependent contact rates
    
    if (Ti <= (Ti + i*h))
        beta_a1 = beta_ia + ((u*beta_ia -beta_ia)/theta)*((Ti+i*h) - Ti);   % Contact rate for asymptomatic cases in the isolated environment
        beta_s1 = beta_is + ((u*beta_is -beta_is)/theta)*((Ti+i*h) - Ti);   % Contact rate for symptomatic cases in the isolated environment

        beta_a2 = beta_ia + ((q1*beta_ia -beta_ia)/theta1)*((Ti+i*h) - Ti); % Contact rate for asymptomatic cases in the non-isolated environment
        beta_s2 = beta_is + ((q1*beta_is -beta_is)/theta1)*((Ti+i*h) - Ti); % Contact rate for symptomatic cases in the non-isolated environment                
    end

    if (Ti + theta <= (Ti + i*h))
        beta_a1 = u*beta_ia;
        beta_s1 = u*beta_is;
    end

    if (Ti + theta1 <= (Ti + i*h))
        beta_a2 = q1*beta_ia;
        beta_s2 = q1*beta_is;
    end

    %% Saving solution
    fprintf(file_with,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.6f\t %1.6f\t %1.6f\t %1.6f\n',(Ti + i*h),s,e,ya,ys,r,c_s,c_e,c_ya,c_ys,c_r,temp_x,beta_a1,beta_s1,beta_a2,beta_s2);

end

fclose(file_with)

%% Prevalence and incidence graphics
NC = 14000000;  % Total population in Wuhan

delimiterIn = ' ';

file_non  = 'dynamics_nonisolation.txt';
file_with = 'dynamics_withisolation.txt';

A_noni = importdata(file_non,delimiterIn);
A_with = importdata(file_with,delimiterIn);

figure
plot(A_noni(:,1),(A_noni(:,4) + A_noni(:,5))*NC,'b')
hold on
plot(A_with(:,1),(A_with(:,4) + A_with(:,5))*NC,'r--')
ylabel('Prevalence')
xlabel('Time')

figure
plot(A_noni(:,1),NC*(A_noni(:,end-1).*A_noni(:,4) + A_noni(:,end).*A_noni(:,5)).*A_noni(:,2),'b')
hold on
plot(A_with(:,1),NC*(A_with(:,end-3).*A_with(:,4) + A_with(:,end-2).*A_with(:,5)).*A_with(:,2),'r--')
ylabel('Incidence')
xlabel('Time')

%% Calculating the maximum value and associated time of the incidence

MaxValue_temp  = max(NC*(A_with(:,end-3).*A_with(:,4) + A_with(:,end-2).*A_with(:,5)).*A_with(:,2));
Index_MaxValue = find((NC*(A_with(:,end-3).*A_with(:,4) + A_with(:,end-2).*A_with(:,5)).*A_with(:,2)) == MaxValue_temp);
Values = [A_with(Index_MaxValue,1)  (NC*(A_with(Index_MaxValue,end-3)*A_with(Index_MaxValue,4) + A_with(Index_MaxValue,end-2)*A_with(Index_MaxValue,5))*A_with(Index_MaxValue,2))]


%% Cumulative incidence vs Hubei data

data_hubei_Jop = [444 549 761 1058 1423 3554 3554 4903 5806 7153 11177 13522 ...
    16678 19665 22112 24953 27100 29631 31728 33366 33366 48206 54406 56249 ...
    58182 59989 61682 62031 62442 62662 64084 64084 64287 64786 65187 65596 ...
    65914 66337 66907 67103 67217 67332 67466 67592 67666 67707 67743 67760 ...
    67773 67781 67786 67790 67794 67798 67799];

time_data = (Ti - 3):1:(size(data_hubei_Jop,2) + Ti - 4);

j1 = 1;
for i=1:size(A_noni,1)
    if A_noni(i,1) == j1
        temp_non(j1) = NC*(A_noni(i,end-1)*A_noni(i,4) + A_noni(i,end)*A_noni(i,5))*A_noni(i,2);
        j1 = j1 + 1;
    end
end

j1 = Ti + 1;
j2 = 1;
for i=1:size(A_with,1)
    if A_with(i,1) == j1
        temp_with_1(j2) = NC*(A_with(i,end-3)*A_with(i,4) + A_with(i,end-2)*A_with(i,5))*A_with(i,2);
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
plot(final_time,final_cum_in,'b')
hold on
plot(time_data,data_hubei_Jop,'r.')
plot(final_time(round(Values(1))),final_cum_in(round(Values(1))),'k*')
plot(time_data(1,8:14),data_hubei_Jop(1,8:14),'k*')
axis([0 100 0 80000])
xlabel('Time')
legend('Cumulative incidence','JHU Data','Location','northwest')

period_time = datetime('now') - initial_time