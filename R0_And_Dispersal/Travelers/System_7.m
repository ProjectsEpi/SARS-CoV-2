% This program solves system 7 of the manuscript "The SARS-CoV-2 epidemic
% outbreak: a review of plausible scenarios of containment and mitigation for
% Mexico" using the 4th order Runge-Kutta method. The outputs are: 
% 1) .txt files with the solution for each instant of time, 2) .txt files with
% estimations of exported cases to some Mexican states, 3) prevalence graph,
% and 4) .txt files of accumulated and daily exported cases from mexico city.
% These last files will be used in the Bar_Graphs.m program to make a bar graph.

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

h11 = 1.25;                 % increase of travelers by holy week
rho_t = 1 - (N_rho/Lim);    % proportion of air infected passengers

%% Basic reproductive number
R0 = sqrt((rho_i*beta_ia + (1 - rho_i)*beta_is)*((s*gamma)/((eta + mu)*(gamma + mu))));

%% Average travelers from Mexico City to some local Mexican States per day
L_mexcan = 6496;
L_mexchi = 1053;
L_mexcju = 974;
L_mexgua = 4535;
L_mexher = 1131;
L_mexmer = 2505;
L_mexmon = 4768;
L_mexoax = 1028;
L_mexpva = 1486;
L_mexsjc = 1341;
L_mextij = 3222;
L_mextux = 1206;
L_mexvil = 1029;

L_mexpue  = 18851;
L_mexgue  = 10799;
L_mexqro  = 6231;
L_mexmic  = 14021;
L_mexhid  = 8728;
L_mextlx  = 3885;
L_mexmor  = 5830;
L_mexguan = 17880;
L_mexjal  = 24027;

%% Open files
file        = fopen('dynamics_travelers.txt','w');

file_cancun = fopen('dynamics_cancun.txt','w');
file_chihua = fopen('dynamics_chihua.txt','w');
file_cjuare = fopen('dynamics_cjuare.txt','w');
file_guadal = fopen('dynamics_guadal.txt','w');
file_hermos = fopen('dynamics_hermos.txt','w');
file_merida = fopen('dynamics_merida.txt','w');
file_monter = fopen('dynamics_monter.txt','w');
file_oaxaca = fopen('dynamics_oaxaca.txt','w');
file_ptoval = fopen('dynamics_ptoval.txt','w');
file_sjcabo = fopen('dynamics_sjcabo.txt','w');
file_tijuan = fopen('dynamics_tijuan.txt','w');
file_tuxlag = fopen('dynamics_tuxlag.txt','w');
file_vilher = fopen('dynamics_vilher.txt','w');

file_puebla = fopen('dynamics_puebla.txt','w');
file_guerre = fopen('dynamics_guerre.txt','w');
file_queret = fopen('dynamics_queret.txt','w');
file_michoa = fopen('dynamics_michoa.txt','w');
file_hidalg = fopen('dynamics_hidalg.txt','w');
file_tlaxca = fopen('dynamics_tlaxca.txt','w');
file_morelo = fopen('dynamics_morelo.txt','w');
file_guanaj = fopen('dynamics_guanaj.txt','w');
file_jalisc = fopen('dynamics_jalisc.txt','w');

%% Exported infected people from Mexico City
lambda_cancun = (L_mexcan)*(ya + ys);
lambda_chihua = (L_mexchi)*(ya + ys);
lambda_cjuare = (L_mexcju)*(ya + ys);
lambda_guadal = (L_mexgua)*(ya + ys);
lambda_hermos = (L_mexher)*(ya + ys);
lambda_merida = (L_mexmer)*(ya + ys);
lambda_monter = (L_mexmon)*(ya + ys);
lambda_oaxaca = (L_mexoax)*(ya + ys);
lambda_ptoval = (L_mexpva)*(ya + ys);
lambda_sjcabo = (L_mexsjc)*(ya + ys);
lambda_tijuan = (L_mextij)*(ya + ys);
lambda_tuxlag = (L_mextux)*(ya + ys);
lambda_vilher = (L_mexvil)*(ya + ys);

lambda_puebla = (L_mexpue)*(ya + ys);
lambda_guerre = (L_mexgue)*(ya + ys);
lambda_queret = (L_mexqro)*(ya + ys);
lambda_michoa = (L_mexmic)*(ya + ys);
lambda_hidalg = (L_mexhid)*(ya + ys);
lambda_tlaxca = (L_mextlx)*(ya + ys);
lambda_morelo = (L_mexmor)*(ya + ys);
lambda_guanaj = (L_mexguan)*(ya + ys);
lambda_jalisc = (L_mexjal)*(ya + ys); 

%%
T   = 200;
tim = 0;        % February 21st
h   = 0.01;     % Step size
N1  = T/h;

%% Saving data
fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',tim,s,e,ya,ys,r,temp_x,N);

fprintf(file_cancun,'%3.6f\t %4.6f\n',tim,lambda_cancun);
fprintf(file_chihua,'%3.6f\t %4.6f\n',tim,lambda_chihua);
fprintf(file_cjuare,'%3.6f\t %4.6f\n',tim,lambda_cjuare);
fprintf(file_guadal,'%3.6f\t %4.6f\n',tim,lambda_guadal);
fprintf(file_hermos,'%3.6f\t %4.6f\n',tim,lambda_hermos);
fprintf(file_merida,'%3.6f\t %4.6f\n',tim,lambda_merida);
fprintf(file_monter,'%3.6f\t %4.6f\n',tim,lambda_monter);
fprintf(file_oaxaca,'%3.6f\t %4.6f\n',tim,lambda_oaxaca);
fprintf(file_ptoval,'%3.6f\t %4.6f\n',tim,lambda_ptoval);
fprintf(file_sjcabo,'%3.6f\t %4.6f\n',tim,lambda_sjcabo);
fprintf(file_tijuan,'%3.6f\t %4.6f\n',tim,lambda_tijuan);
fprintf(file_tuxlag,'%3.6f\t %4.6f\n',tim,lambda_tuxlag);
fprintf(file_vilher,'%3.6f\t %4.6f\n',tim,lambda_vilher);

fprintf(file_puebla,'%3.6f\t %4.6f\n',tim,lambda_puebla);
fprintf(file_guerre,'%3.6f\t %4.6f\n',tim,lambda_guerre);
fprintf(file_queret,'%3.6f\t %4.6f\n',tim,lambda_queret);
fprintf(file_michoa,'%3.6f\t %4.6f\n',tim,lambda_michoa);
fprintf(file_hidalg,'%3.6f\t %4.6f\n',tim,lambda_hidalg);
fprintf(file_tlaxca,'%3.6f\t %4.6f\n',tim,lambda_tlaxca);
fprintf(file_morelo,'%3.6f\t %4.6f\n',tim,lambda_morelo);
fprintf(file_guanaj,'%3.6f\t %4.6f\n',tim,lambda_guanaj);
fprintf(file_jalisc,'%3.6f\t %4.6f\n',tim,lambda_jalisc);

%% Solving EDO System with 4th order Runge-Kutta method

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
    
    %%
    
    if (43 <= (tim + h)) && ((tim + h) <= 51)
       
        % Average travelers from Mexico City to some local Mexican States
        % per day (Holy Week)
        L_mexcan = 6496*h11;
        L_mexchi = 1053*h11;
        L_mexcju = 974*h11;
        L_mexgua = 4535*h11;
        L_mexher = 1131*h11;
        L_mexmer = 2505*h11;
        L_mexmon = 4768*h11;
        L_mexoax = 1028*h11;
        L_mexpva = 1486*h11;
        L_mexsjc = 1341*h11;
        L_mextij = 3222*h11;
        L_mextux = 1206*h11;
        L_mexvil = 1029*h11;         
        
        L_mexpue  = 18851*h11;
        L_mexgue  = 10799*h11;
        L_mexqro  = 6231*h11;
        L_mexmic  = 14021*h11;
        L_mexhid  = 8728*h11;
        L_mextlx  = 3885*h11;
        L_mexmor  = 5830*h11;
        L_mexguan = 17880*h11;
        L_mexjal  = 24027*h11;
        
    else
       
        % Average travelers from Mexico City to some local Mexican States
        % per day (except Holy Week)
        L_mexcan = 6496;
        L_mexchi = 1053;
        L_mexcju = 974;
        L_mexgua = 4535;
        L_mexher = 1131;
        L_mexmer = 2505;
        L_mexmon = 4768;
        L_mexoax = 1028;
        L_mexpva = 1486;
        L_mexsjc = 1341;
        L_mextij = 3222;
        L_mextux = 1206;
        L_mexvil = 1029;
        
        L_mexpue  = 18851;
        L_mexgue  = 10799;
        L_mexqro  = 6231;
        L_mexmic  = 14021;
        L_mexhid  = 8728;
        L_mextlx  = 3885;
        L_mexmor  = 5830;
        L_mexguan = 17880;
        L_mexjal  = 24027;        
    end
    
    %% Exported infected people from Mexico City
    lambda_cancun = (L_mexcan)*(ya + ys);
    lambda_chihua = (L_mexchi)*(ya + ys);
    lambda_cjuare = (L_mexcju)*(ya + ys);
    lambda_guadal = (L_mexgua)*(ya + ys);
    lambda_hermos = (L_mexher)*(ya + ys);
    lambda_merida = (L_mexmer)*(ya + ys);
    lambda_monter = (L_mexmon)*(ya + ys);
    lambda_oaxaca = (L_mexoax)*(ya + ys);
    lambda_ptoval = (L_mexpva)*(ya + ys);
    lambda_sjcabo = (L_mexsjc)*(ya + ys);
    lambda_tijuan = (L_mextij)*(ya + ys);
    lambda_tuxlag = (L_mextux)*(ya + ys);
    lambda_vilher = (L_mexvil)*(ya + ys);

    lambda_puebla = (L_mexpue)*(ya + ys);
    lambda_guerre = (L_mexgue)*(ya + ys);
    lambda_queret = (L_mexqro)*(ya + ys);
    lambda_michoa = (L_mexmic)*(ya + ys);
    lambda_hidalg = (L_mexhid)*(ya + ys);
    lambda_tlaxca = (L_mextlx)*(ya + ys);
    lambda_morelo = (L_mexmor)*(ya + ys);
    lambda_guanaj = (L_mexguan)*(ya + ys);
    lambda_jalisc = (L_mexjal)*(ya + ys);    

    %% Saving data
    fprintf(file,'%3.6f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %1.9f\t %9.9f\n',tim + h,s,e,ya,ys,r,temp_x,N);
    
    fprintf(file_cancun,'%3.6f\t %4.6f\n',tim + h,lambda_cancun);
    fprintf(file_chihua,'%3.6f\t %4.6f\n',tim + h,lambda_chihua);
    fprintf(file_cjuare,'%3.6f\t %4.6f\n',tim + h,lambda_cjuare);
    fprintf(file_guadal,'%3.6f\t %4.6f\n',tim + h,lambda_guadal);
    fprintf(file_hermos,'%3.6f\t %4.6f\n',tim + h,lambda_hermos);
    fprintf(file_merida,'%3.6f\t %4.6f\n',tim + h,lambda_merida);
    fprintf(file_monter,'%3.6f\t %4.6f\n',tim + h,lambda_monter);
    fprintf(file_oaxaca,'%3.6f\t %4.6f\n',tim + h,lambda_oaxaca);
    fprintf(file_ptoval,'%3.6f\t %4.6f\n',tim + h,lambda_ptoval);
    fprintf(file_sjcabo,'%3.6f\t %4.6f\n',tim + h,lambda_sjcabo);
    fprintf(file_tijuan,'%3.6f\t %4.6f\n',tim + h,lambda_tijuan);
    fprintf(file_tuxlag,'%3.6f\t %4.6f\n',tim + h,lambda_tuxlag);
    fprintf(file_vilher,'%3.6f\t %4.6f\n',tim + h,lambda_vilher);
    
    fprintf(file_puebla,'%3.6f\t %4.6f\n',tim + h,lambda_puebla);
    fprintf(file_guerre,'%3.6f\t %4.6f\n',tim + h,lambda_guerre);
    fprintf(file_queret,'%3.6f\t %4.6f\n',tim + h,lambda_queret);
    fprintf(file_michoa,'%3.6f\t %4.6f\n',tim + h,lambda_michoa);
    fprintf(file_hidalg,'%3.6f\t %4.6f\n',tim + h,lambda_hidalg);
    fprintf(file_tlaxca,'%3.6f\t %4.6f\n',tim + h,lambda_tlaxca);
    fprintf(file_morelo,'%3.6f\t %4.6f\n',tim + h,lambda_morelo);
    fprintf(file_guanaj,'%3.6f\t %4.6f\n',tim + h,lambda_guanaj);
    fprintf(file_jalisc,'%3.6f\t %4.6f\n',tim + h,lambda_jalisc);  

end

%% Close files
fclose(file);

fclose(file_cancun);
fclose(file_chihua);
fclose(file_cjuare);
fclose(file_guadal);
fclose(file_hermos);
fclose(file_merida);
fclose(file_monter);
fclose(file_oaxaca);
fclose(file_ptoval);
fclose(file_sjcabo);
fclose(file_tijuan);
fclose(file_tuxlag);
fclose(file_vilher);

fclose(file_puebla);
fclose(file_guerre);
fclose(file_queret);
fclose(file_michoa);
fclose(file_hidalg);
fclose(file_tlaxca);
fclose(file_morelo);
fclose(file_guanaj);
fclose(file_jalisc);

%% Prevalence Graphic
delimiterIn = ' ';
file_trav   = 'dynamics_travelers.txt';
A_trav      = importdata(file_trav,delimiterIn);

figure
plot(A_trav(:,1),(A_trav(:,4) + A_trav(:,5)),'b')
ylabel('Prevalence')
xlabel('Time')

%% Exported cases estimation accumulated and per day
file_cancun = 'dynamics_cancun.txt';
file_chihua = 'dynamics_chihua.txt';
file_cjuare = 'dynamics_cjuare.txt';
file_guadal = 'dynamics_guadal.txt';
file_hermos = 'dynamics_hermos.txt';
file_merida = 'dynamics_merida.txt';
file_monter = 'dynamics_monter.txt';
file_oaxaca = 'dynamics_oaxaca.txt';
file_ptoval = 'dynamics_ptoval.txt';
file_sjcabo = 'dynamics_sjcabo.txt';
file_tijuan = 'dynamics_tijuan.txt';
file_tuxlag = 'dynamics_tuxlag.txt';
file_vilher = 'dynamics_vilher.txt';

file_puebla = 'dynamics_puebla.txt';
file_guerre = 'dynamics_guerre.txt';
file_queret = 'dynamics_queret.txt';
file_michoa = 'dynamics_michoa.txt';
file_hidalg = 'dynamics_hidalg.txt';
file_tlaxca = 'dynamics_tlaxca.txt';
file_morelo = 'dynamics_morelo.txt';
file_guanaj = 'dynamics_guanaj.txt';
file_jalisc = 'dynamics_jalisc.txt';

A_cancun = importdata(file_cancun,delimiterIn);
A_chihua = importdata(file_chihua,delimiterIn);
A_cjuare = importdata(file_cjuare,delimiterIn);
A_guadal = importdata(file_guadal,delimiterIn);
A_hermos = importdata(file_hermos,delimiterIn);
A_merida = importdata(file_merida,delimiterIn);
A_monter = importdata(file_monter,delimiterIn);
A_oaxaca = importdata(file_oaxaca,delimiterIn);
A_ptoval = importdata(file_ptoval,delimiterIn);
A_sjcabo = importdata(file_sjcabo,delimiterIn);
A_tijuan = importdata(file_tijuan,delimiterIn);
A_tuxlag = importdata(file_tuxlag,delimiterIn);
A_vilher = importdata(file_vilher,delimiterIn);

A_puebla = importdata(file_puebla,delimiterIn);
A_guerre = importdata(file_guerre,delimiterIn);
A_queret = importdata(file_queret,delimiterIn);
A_michoa = importdata(file_michoa,delimiterIn);
A_hidalg = importdata(file_hidalg,delimiterIn);
A_tlaxca = importdata(file_tlaxca,delimiterIn);
A_morelo = importdata(file_morelo,delimiterIn);
A_guanaj = importdata(file_guanaj,delimiterIn);
A_jalisc = importdata(file_jalisc,delimiterIn);

file_puebla = 'dynamics_puebla.txt';
file_guerre = 'dynamics_guerre.txt';
file_queret = 'dynamics_queret.txt';
file_michoa = 'dynamics_michoa.txt';
file_hidalg = 'dynamics_hidalg.txt';
file_tlaxca = 'dynamics_tlaxca.txt';
file_morelo = 'dynamics_morelo.txt';
file_guanaj = 'dynamics_guanaj.txt';
file_jalisc = 'dynamics_jalisc.txt';

file1 = fopen('cumulative_air.txt','w');
file2 = fopen('cumulative_road.txt','w');
file3 = fopen('PerDay_air.txt','w');
file4 = fopen('PerDay_road.txt','w');

j_temp = 1;
k_temp = 1;

sum_cancun = 0;
sum_chihua = 0;
sum_cjuare = 0;
sum_guadal = 0;
sum_hermos = 0;
sum_merida = 0;
sum_monter = 0;
sum_oaxaca = 0;
sum_ptoval = 0;
sum_sjcabo = 0;
sum_tijuan = 0;
sum_tuxlag = 0;
sum_vilher = 0;

sum_puebla = 0;
sum_guerre = 0;
sum_queret = 0;
sum_michoa = 0;
sum_hidalg = 0;
sum_tlaxca = 0;
sum_morelo = 0;
sum_guanaj = 0;
sum_jalisc = 0;

for i = 1:size(A_cancun,1)
    if A_cancun(i,1) == k_temp
        
        Q1_cancun = h*trapz(A_cancun(j_temp:i,2));
        Q1_chihua = h*trapz(A_chihua(j_temp:i,2));
        Q1_cjuare = h*trapz(A_cjuare(j_temp:i,2));
        Q1_guadal = h*trapz(A_guadal(j_temp:i,2));
        Q1_hermos = h*trapz(A_hermos(j_temp:i,2));
        Q1_merida = h*trapz(A_merida(j_temp:i,2));
        Q1_monter = h*trapz(A_monter(j_temp:i,2));
        Q1_oaxaca = h*trapz(A_oaxaca(j_temp:i,2));
        Q1_ptoval = h*trapz(A_ptoval(j_temp:i,2));
        Q1_sjcabo = h*trapz(A_sjcabo(j_temp:i,2));
        Q1_tijuan = h*trapz(A_tijuan(j_temp:i,2));
        Q1_tuxlag = h*trapz(A_tuxlag(j_temp:i,2));
        Q1_vilher = h*trapz(A_vilher(j_temp:i,2));
        
        Q1_puebla = h*trapz(A_puebla(j_temp:i,2));
        Q1_guerre = h*trapz(A_guerre(j_temp:i,2));
        Q1_queret = h*trapz(A_queret(j_temp:i,2));
        Q1_michoa = h*trapz(A_michoa(j_temp:i,2));
        Q1_hidalg = h*trapz(A_hidalg(j_temp:i,2));
        Q1_tlaxca = h*trapz(A_tlaxca(j_temp:i,2)); 
        Q1_morelo = h*trapz(A_morelo(j_temp:i,2));
        Q1_guanaj = h*trapz(A_guanaj(j_temp:i,2));
        Q1_jalisc = h*trapz(A_jalisc(j_temp:i,2));

        %% Cumulative Exported Cases
        sum_cancun = Q1_cancun + sum_cancun;
        sum_chihua = Q1_chihua + sum_chihua;
        sum_cjuare = Q1_cjuare + sum_cjuare;
        sum_guadal = Q1_guadal + sum_guadal;
        sum_hermos = Q1_hermos + sum_hermos;
        sum_merida = Q1_merida + sum_merida;
        sum_monter = Q1_monter + sum_monter;
        sum_oaxaca = Q1_oaxaca + sum_oaxaca;
        sum_ptoval = Q1_ptoval + sum_ptoval;
        sum_sjcabo = Q1_sjcabo + sum_sjcabo;
        sum_tijuan = Q1_tijuan + sum_tijuan;
        sum_tuxlag = Q1_tuxlag + sum_tuxlag;
        sum_vilher = Q1_vilher + sum_vilher;

        sum_puebla = Q1_puebla + sum_puebla;
        sum_guerre = Q1_guerre + sum_guerre;
        sum_queret = Q1_queret + sum_queret;
        sum_michoa = Q1_michoa + sum_michoa;
        sum_hidalg = Q1_hidalg + sum_hidalg;
        sum_tlaxca = Q1_tlaxca + sum_tlaxca;
        sum_morelo = Q1_morelo + sum_morelo;
        sum_guanaj = Q1_guanaj + sum_guanaj; 
        sum_jalisc = Q1_jalisc + sum_jalisc;        
        
        
        fprintf(file1,'%3.0f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\n',...
            k_temp,sum_cancun,sum_chihua,sum_cjuare,sum_guadal,sum_hermos,sum_merida,sum_monter,sum_oaxaca,sum_ptoval,sum_sjcabo,sum_tijuan,...
            sum_tuxlag,sum_vilher);      
        
        fprintf(file2,'%3.0f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\n',k_temp,sum_puebla,sum_guerre,sum_queret,...
            sum_michoa,sum_hidalg,sum_tlaxca,sum_morelo,sum_guanaj,sum_jalisc); 
        
        fprintf(file3,'%3.0f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\n',...
            k_temp,Q1_cancun,Q1_chihua,Q1_cjuare,Q1_guadal,Q1_hermos,Q1_merida,Q1_monter,Q1_oaxaca,Q1_ptoval,Q1_sjcabo,Q1_tijuan,...
            Q1_tuxlag,Q1_vilher);      
        
        fprintf(file4,'%3.0f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\t %5.6f\n',k_temp,Q1_puebla,Q1_guerre,Q1_queret,...
            Q1_michoa,Q1_hidalg,Q1_tlaxca,Q1_morelo,Q1_guanaj,Q1_jalisc);         
        
        j_temp = i + 1;
        k_temp = k_temp + 1;
    end
end

fclose(file1)
fclose(file2)
fclose(file3)
fclose(file4)

period_time = datetime('now') - initial_time