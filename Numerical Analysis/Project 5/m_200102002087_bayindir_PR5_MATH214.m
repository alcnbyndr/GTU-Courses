%% Alican Bayındır 200102002087
% MATH 214 - Project 5
% 05.01.2021
close all; clear all; clc;
format short;
load pr5data.dat;

kilometer_values = pr5data(:, 1)'./1000;
voltage_values = pr5data(:, 2)';

x = 0; x2 = 0; y = 0; xy = 0; elsize = 11;

for i = 1:elsize
    x = x + (kilometer_values(1,i));
    y = y + (voltage_values(1,i));  
    xy = xy + (kilometer_values(1,i) * voltage_values(1,i));
    x2 = x2 + (kilometer_values(1,i) * kilometer_values(1,i));
end

a0 = ((x2 * y) - (xy * x)) / ((elsize * x2) - ((x) * (x)));
a1 = ((elsize * xy) - (x * y)) / ((elsize * x2) - ((x) * (x)));

yi = [];
for i = 1:elsize
    yi = [yi (a1 * kilometer_values(1,i) + a0)];
end

x3 = 0;
x4 = 0;
yx2 = 0;
for i = 1:elsize
    x3 = x3 + (kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));
    x4 = x4 + (kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));
    yx2 = yx2 + (voltage_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));
end

A = [elsize x x2; x x2 x3; x2 x3 x4];
S = [y; xy; yx2];
B = inv(A) * S;

Y = [];
for i = 1:elsize
    Y = [Y ((B(3,1) * kilometer_values(1,i) * kilometer_values(1,i)) + (B(2,1) * kilometer_values(1,i) + (B(1,1))))];
end

plot(kilometer_values, voltage_values, '-*', 'LineWidth', 2);
hold on;

plot(kilometer_values, yi, '-*', 'LineWidth', 2);
hold on;

plot(kilometer_values, Y,'-*', 'LineWidth', 2);
xlabel('Kilometer'); ylabel('Voltages'); grid on;
title('Voltage-Kilometer graph');
hold on;

x5 = 0;
for i = 1:elsize
    x5 = x5 + (kilometer_values(1,i) * kilometer_values(1,i) * ...
    kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));
end

x6=0;
for i = 1:elsize
    x6 = x6 + (kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i) ...
    * kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));
end

yx3 = 0;
for i = 1:elsize
    yx3 = yx3 + (voltage_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i) * kilometer_values(1,i));   
end

A2 = [elsize x x2 x3; x x2 x3 x4; x2 x3 x4 x5; x3 x4 x5 x6];
S2 = [y; xy; yx2; yx3];
B2 = inv(A2) * S2;

Y2  = [];
for i = 1:elsize 
    Y2 = [Y2 ((B2(4,1) * kilometer_values(1,i) * kilometer_values(1,i) * ...
    kilometer_values(1,i)) + (B2(3,1) * kilometer_values(1,i) * ...
    kilometer_values(1,i)) + (B2(2,1) * kilometer_values(1,i)+(B2(1,1))))];
end
plot(kilometer_values, Y2,'-*', 'LineWidth', 2);
legend('Real Value', 'Linear', 'Degree 2', 'Degree 3');

error_degree_linear = []; error_degree_two = []; error_degree_three = [];
for i = 1:elsize
    error_degree_linear = [error_degree_linear (voltage_values(1,i)-yi(1,i))];
    error_degree_two = [error_degree_two (voltage_values(1,i)-Y(1,i))];
    error_degree_three = [error_degree_three (voltage_values(1,i)-Y2(1,i))];
end