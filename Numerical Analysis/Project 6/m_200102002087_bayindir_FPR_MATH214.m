%% Alican Bayındır 200102002087
% MATH 214 - Final Project
% 22.01.2021
clc; clear; close all;
format long;

% Load the data
load fprdata.dat

% Edit the data to be used easily in future operations
voltage_values = fprdata(:, 1)';
current_values = fprdata(:, 2)';

% Constant values given in the project guideline
VS = 2; L = 0.98; R = 14.2;
DELTA_T1 = 0.025; DELTA_T2 = 0.0025; 

% It is required to define the part below since we use it in Euler's Method
% to make for loop run according to STEP_COUNT numbers.
STEP_COUNT = 0:DELTA_T1:0.6;
STEP_COUNT2 = 0:DELTA_T2:0.6;
elsize = length(fprdata);

% The part below iss defined for Euler's Method.
derivative_of_current1(1) = VS/L; derivative_of_current2(1) = VS/L;
current_data1(1) = 0; current_data2(1) = 0;
step_one(1)  =  0; step_two(1) = 0;

VD_one(1) = 0; VD_two(1) = 0;
VR_one(1) = 0; VR_two(1) = 0;
VL_one(1) = VS; VL_two(1) = VS;

% The part below takes the log of the current values that are given
% in the fprdata.dat file.
for i = 1:elsize
    log_of_current(i) = log10(current_values(i));
end

% The following values is needed to find a, b and lnb values.
ln_current = 0; sqrt_current = 0; ln_currents = 0; sum_of_voltages = 0;
for i = 1:elsize
    sum_of_voltages  =  sum_of_voltages + voltage_values(i);
    ln_current  =  ln_current + log(current_values(i));
    sqrt_current  =  sqrt_current + voltage_values(i).^2;
    ln_currents  =  ln_currents + voltage_values(i) * log(current_values(i));
end

% These values will be used to calculate The Exponential form of Least 
% Square Approximation.
a  =  (5 * ln_currents - sum_of_voltages * ln_current) / (5 * sqrt_current - sum_of_voltages.^2);
lnb  =  (sqrt_current * ln_current - ln_currents * sum_of_voltages) / (5 * sqrt_current - sum_of_voltages.^2);
b  =  exp(lnb);

% The part below is the definition of the equation 1.3 in the project
% report in MATLAB platform.
for i = 1:elsize
    current(i) =  b * exp(voltage_values(i) * a);
    function_log(i) =  log10(current(i));
end

% Euler's method's definition when the DELTA_T1 equals 0.025 ms
for i = 2:length(STEP_COUNT)
    current_data1(i) = current_data1(i-1) + DELTA_T1 * derivative_of_current1(i-1);
    VD_one(i) = (log(current_data1(i)) - log(b)) / a;
    derivative_of_current1(i) = (VS - current_data1(i) * R - VD_one(i)) / L;
    step_one(i) = step_one(i-1) + DELTA_T1;
    VR_one(i) =  current_data1(i) * R;
    VL_one(i) = derivative_of_current1(i) * L;
end

% Euler's method's definition when the DELTA_T2 equals 0.0025 ms
for i = 2:length(STEP_COUNT2)
    current_data2(i) = current_data2(i-1) + DELTA_T2 * derivative_of_current2(i-1);
    VD_two(i) = (log(current_data2(i)) - log(b) ) / a;
    derivative_of_current2(i) = (VS - current_data2(i)*R - VD_two(i)) / L;
    step_two(i) = step_two(i-1) + DELTA_T2;
    VR_two(i) = current_data2(i) * R;
    VL_two(i) = derivative_of_current2(i) * L;
end

% Plotting the Current-Voltage graph's values that are given in the fprdata.dat file 
% and the current values that is calculated by using equation (1.3) in the
% report.
figure(1);
plot(voltage_values, current_values, "rd", voltage_values, current, "-*", 'LineWidth', 2);
title('Voltage - Current Graph'); xlabel('Voltage [V]'); ylabel('Current [A]');
legend('Data from fprdata.dat file','Approximated values', 'Location', 'northwest');
grid on;

% Plotting the Log of the current and log of the function with respect to voltage values.
figure(2);
plot(voltage_values, log_of_current, "rd", voltage_values, function_log, "-*", 'LineWidth', 2);
title('Voltage-Current Graph in Logarithmic scale'); xlabel('Voltage [V]'); ylabel('Current [A]');
legend('Data from fprdata.dat file','Approximated values', 'Location', 'northwest');
grid on;

% The subplots below is the plotting part of the Euler's method's values.
figure(3);
subplot(2,2,1);
plot(step_one, current_data1, "-x", step_two, current_data2, "-x");
title('Current i[t] Graph'); xlabel('Time'); ylabel('Current [A]');
legend('Step size 0.025','Step size 0.0025', 'Location', 'southeast');
grid on;

subplot(2,2,2);
plot(step_one, VD_one, "-x", step_two, VD_two, "-x");
title('Diode Voltage [VD] Graph'); xlabel('Time'); ylabel('Voltage [V]');
legend('Step size 0.025','Step size 0.0025', 'Location', 'southeast');
grid on;

subplot(2,2,3);
plot(step_one, VR_one, "-x", step_two, VR_two, "-x");
title('Resistance Voltage [VR] Graph'); xlabel('Time'); ylabel('Voltage [V]');
legend('Step size 0.025','Step size 0.0025', 'Location', 'southeast');
grid on;

subplot(2,2,4);
plot(step_one, VL_one, "-x", step_two, VL_two, "-x");
title('Inductance Voltage [VL] Graph'); xlabel('Time'); ylabel('Voltage [V]');
legend('Step size 0.025','Step size 0.0025');
grid on;