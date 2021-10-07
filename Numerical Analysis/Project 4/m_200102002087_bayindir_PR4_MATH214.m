%% Alican Bayındır 200102002087
% MATH 214 - Project 4
% 23.12.2020
close all; clear all; clc; 

STEP_SIZE_1 = 0.05;
STEP_SIZE_2 = 0.025;
L = 0.98; R = 14.2; Vs = 12;
INITIAL_CURRENT = 0.1;

% To draw plot we need to specify the time values
time_half = 0:STEP_SIZE_1:0.6;
time_quarter = 0:STEP_SIZE_2:0.6;

% Some initial values that should be specified before starting methods
current_euler_half(1) = INITIAL_CURRENT;
current_euler_quarter(1) = INITIAL_CURRENT;

current_modieuler_half(1) = INITIAL_CURRENT;
current_modieuler_quarter(1) = INITIAL_CURRENT;

current_midpoint_half(1) = INITIAL_CURRENT;
current_midpoint_quarter(1) = INITIAL_CURRENT;

current_runge_kutta_half(1) = INITIAL_CURRENT;
current_runge_kutta_quarter(1) = INITIAL_CURRENT;

current_analytical_half(1) = INITIAL_CURRENT;
current_analytical_quarter(1) = INITIAL_CURRENT;

F = @(y) (Vs - R * y) / L;
current_equation = @(t) ((Vs*(1-exp((-R*t)/L))) / R) + INITIAL_CURRENT;

% Euler's method
% Step size = 0.05;
for k = 1:length(time_half)-1
    current_euler_half(k+1) = current_euler_half(k) + STEP_SIZE_1 * F(current_euler_half(k));
end

% Step size = 0.025;
for k = 1:length(time_quarter)-1
    current_euler_quarter(k+1) = current_euler_quarter(k) + STEP_SIZE_2 * F(current_euler_quarter(k));
end

% Modified Euler's Method
for k = 1:length(time_half)-1
    current_modieuler_half(k+1) = current_modieuler_half(k) + (STEP_SIZE_1 / 2) * (F(current_modieuler_half(k)) + STEP_SIZE_1 * F(current_modieuler_half(k)));
end

% Modified Euler's Method
for k = 1:length(time_quarter)-1
    current_modieuler_quarter(k+1) = current_modieuler_quarter(k) + (STEP_SIZE_2 / 2) * (F(current_modieuler_quarter(k)) + STEP_SIZE_2 * F(current_modieuler_quarter(k)));
end

% Midpoint Method
for k = 1:length(time_half)-1
    current_midpoint_half(k+1)= current_midpoint_half(k) + STEP_SIZE_1 * (F(current_midpoint_half(k)) + (STEP_SIZE_1 / 2) * F(current_midpoint_half(k)));
end

% Midpoint Method
for k = 1:length(time_quarter)-1
    current_midpoint_quarter(k+1)= current_midpoint_quarter(k) + STEP_SIZE_2 * (F(current_midpoint_quarter(k)) + (STEP_SIZE_2 / 2) * F(current_midpoint_quarter(k)));
end

% Runge Kutta Method Fourth Order
% Step size = 0.05;
for k = 1:length(time_half)-1
    runge_kutta_half1 = STEP_SIZE_1 * F(current_runge_kutta_half(k));
    runge_kutta_half2 = STEP_SIZE_1 * F(current_runge_kutta_half(k) + runge_kutta_half1 / 2);
    runge_kutta_half3 = STEP_SIZE_1 * F(current_runge_kutta_half(k) + runge_kutta_half2 / 2);
    runge_kutta_half4 = STEP_SIZE_1 * F(current_runge_kutta_half(k) + runge_kutta_half3);
    current_runge_kutta_half(k + 1) = current_runge_kutta_half(k) + (runge_kutta_half1 + 2 * runge_kutta_half2 + 2 * runge_kutta_half3 + runge_kutta_half4) / 6;
end

% Step size = 0.025
for k = 1:length(time_quarter)-1
    runge_kutta_quarter1 = STEP_SIZE_2 * F(current_runge_kutta_quarter(k));
    runge_kutta_quarter2 = STEP_SIZE_2 * F(current_runge_kutta_quarter(k) + runge_kutta_quarter1 / 2);
    runge_kutta_quarter3 = STEP_SIZE_2 * F(current_runge_kutta_quarter(k) + runge_kutta_quarter2 / 2);
    runge_kutta_quarter4 = STEP_SIZE_2 * F(current_runge_kutta_quarter(k) + runge_kutta_quarter3);
    current_runge_kutta_quarter(k + 1) = current_runge_kutta_quarter(k) + (runge_kutta_quarter1 + 2 * runge_kutta_quarter2 + 2 * runge_kutta_quarter3 + runge_kutta_quarter4) / 6;
end

% Error Analysis
for i = 1:length(time_half)-1
    current_analytical_half(i+1) = current_equation(STEP_SIZE_1*i);
end
     
for i = 1:length(time_quarter)-1
    current_analytical_quarter(i+1) = current_equation(STEP_SIZE_2*i);
end

% Plots of the methods
figure(1);
plot(time_half,current_analytical_half, '-*', time_half, current_euler_half, '-*', time_half, current_modieuler_half, '-*', time_half, current_midpoint_half, '-*', time_half, current_runge_kutta_half, '-*', 'LineWidth', 2);
xlabel('Time'); ylabel('Current values'); grid on;
title('The graph of current values of each method when delta t = 0.05');
legend('Real','Euler', 'Modified', 'Midpoint', 'Runge Kutta', 'Location', 'southeast');

figure(2);
plot(time_quarter, current_analytical_quarter, '-*', time_quarter, current_euler_quarter, '-*', time_quarter, current_modieuler_quarter, '-*', time_quarter, current_midpoint_quarter, '-*', time_quarter, current_runge_kutta_quarter, '-*', 'LineWidth', 2);
xlabel('Time'); ylabel('Current values'); grid on;
title('The graph of current values of each method when delta t = 0.025');
legend('Real','Euler', 'Modified', 'Midpoint', 'Runge Kutta', 'Location', 'southeast');

figure(3);
subplot(2,2,1);
plot(time_half,current_analytical_half, '-*', time_half, current_euler_half, '-*', time_quarter, current_euler_quarter, '-*');
grid on;
xlabel('Time'); ylabel('Results of Euler Method'); grid on;
title('Euler Method Graph');
legend('Analytical values', 'Delta T = 0.05s', 'Delta T = 0.025s', 'Location', 'southeast');

subplot(2,2,2);
plot(time_half,current_analytical_half, '-*', time_half, current_modieuler_half, '-*', time_quarter, current_modieuler_quarter, '-*');
grid on;
xlabel('Time'); ylabel('Results of Modified Euler Method'); grid on;
title('Modified Euler Method Graph');
legend('Analytical values', 'Delta T = 0.05s', 'Delta T = 0.025s', 'Location', 'southeast');

subplot(2,2,3);
plot(time_half,current_analytical_half, '-*', time_half, current_midpoint_half, '-*', time_quarter, current_midpoint_quarter, '-*');
grid on;
xlabel('Time'); ylabel('Results of Midpoint Method'); grid on;
title('Midpoint Method Graph');
legend('Analytical values', 'Delta T = 0.05s', 'Delta T = 0.025s', 'Location', 'southeast');

subplot(2,2,4);
plot(time_half,current_analytical_half, '-*', time_half, current_runge_kutta_half, '-*', time_quarter, current_runge_kutta_quarter, '-*');
xlabel('Time'); ylabel('Results of Runge-Kutta Method (4th order)'); grid on;
title('Runge-Kutta Method Graph');
legend('Analytical values', 'Delta T = 0.05s', 'Delta T = 0.025s', 'Location', 'southeast');

error_euler_half = abs(current_analytical_half - current_euler_half);
error_modieuler_half=abs(current_analytical_half - current_modieuler_half);
error_midpoint_half=abs(current_analytical_half - current_midpoint_half);
error_rungekutta_half=abs(current_analytical_half - current_runge_kutta_half);

error_euler_quarter=abs(current_analytical_quarter - current_euler_quarter);
error_modieuler_quarter=abs(current_analytical_quarter - current_modieuler_quarter);
error_midpoint_quarter=abs(current_analytical_quarter- current_midpoint_quarter);
error_rungekutta_quarter=abs(current_analytical_quarter - current_runge_kutta_quarter);

figure(4);
subplot(1,2,1);
plot(time_half, error_euler_half, '-*', time_half, error_modieuler_half, '-*', time_half, error_midpoint_half, '-*', time_half, error_rungekutta_half, '-*');
xlabel('Error rate'); ylabel('Time'); grid on;
title('Error rates of all methods when Delta t = 0.05s');
legend('Euler Error', 'Modified Euler Error', 'Midpoint Error', 'Runge-Kutta Error', 'Location', 'southeast');

subplot(1,2,2);
plot(time_quarter, error_euler_quarter, '-*', time_quarter, error_modieuler_quarter, '-*', time_quarter, error_midpoint_quarter, '-*', time_quarter, error_rungekutta_quarter, '-*');
xlabel('Error rate'); ylabel('Time'); grid on;
title('Error rates of all methods when Delta t = 0.025s');
legend('Euler Error', 'Modified Euler Error', 'Midpoint Error', 'Runge-Kutta Error', 'Location', 'southeast');