% This program has as outputs bar graphs for a specified date. These can be
% cumulative data or for a specific day.

clc
clear all

delimiterIn = ' ';

% file_1 = 'PerDay_air_1.txt';
% file_2 = 'PerDay_air_2.txt';
% file_3 = 'PerDay_air_3.txt';
% file_4 = 'PerDay_road_1.txt';
% file_5 = 'PerDay_road_2.txt';
% file_6 = 'PerDay_road_3.txt';

file_1 = 'cumulative_air_1.txt';
file_2 = 'cumulative_air_2.txt';
file_3 = 'cumulative_air_3.txt';
file_4 = 'cumulative_road_1.txt';
file_5 = 'cumulative_road_2.txt';
file_6 = 'cumulative_road_3.txt';

A1 = importdata(file_1,delimiterIn);
A2 = importdata(file_2,delimiterIn);
A3 = importdata(file_3,delimiterIn);
A4 = importdata(file_4,delimiterIn);
A5 = importdata(file_5,delimiterIn);
A6 = importdata(file_6,delimiterIn);

%% Date
date_1 = 60;        % Date to calculate bar graphs

n = size(A1,1);

j1 = 1;

for i = 1:n
    %%
    if A1(i,1) == date_1
        B1(j1,:) = A1(i,:);
        B2(j1,:) = A4(i,:);
        j1 = j1 + 1;
    end    
    
    if A2(i,1) == date_1
        B1(j1,:) = A2(i,:);
        B2(j1,:) = A5(i,:);
        j1 = j1 + 1;
    end
    
    if A3(i,1) == date_1
        B1(j1,:) = A3(i,:);
        B2(j1,:) = A6(i,:);
        j1 = j1 + 1;
    end    
    
end

C1_temp = B1;
C1_temp(:,1) = [];
C1 = C1_temp';
n1 = size(C1,1);

C2_temp = B2;
C2_temp(:,1) = [];
C2 = C2_temp';
n2 = size(C2,1);

j1 = 1;

for i = 1:n1
    D1(i,:) = C1(j1,:);
    j1 = j1 + 1;
end

j2 = 1;
for i = 1:n2
    D2(i,:) = C2(j2,:);
    j2 = j2 + 1;
end

c = categorical({'Cancún','Chihuahua','Ciudad Juarez','Guadalajara','Hermosillo','Mérida','Monterrey','Oaxaca','Pto Vallarta','San José del Cabo','Tijuana','Tuxtla Gutierrez','Villahermosa'});
d = categorical({'Puebla','Guerrero','Querétaro','Michoacán','Hidalgo','Tlaxcala','Morelos','Guanajuato','Jalisco'});

figure
bar(c,D1)

figure
bar(d,D2)