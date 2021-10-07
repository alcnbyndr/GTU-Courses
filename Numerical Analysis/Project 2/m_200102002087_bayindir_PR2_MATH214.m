clear all;
close all;
clc;

% to decrease truncation error I used format long.
format long

L= 0.98; % henry
R= 14.2; % ohm

% Delta values it shows how often we measured the currents.
delta_t1 = 0.075; 
delta_t2 = 0.050;
delta_t3 = 0.025;
delta_t4 = 0.010;

% Reading the datas 
datafiles = {'current1.dat' 'current2.dat' 'current3.dat' 'current4.dat'};

for i = 1:numel(datafiles)
    load( datafiles{i} )
end

% to prevent matrix overflow we need to decrease number of total elements
% in the matrix minus one. 
% the format of the variables are type_file_method so if you see
% derivat_cur2_forward it means it is derivative of the current 2 function
% to apply forward difference method.
% Then, calculated the voltage in the circuit that can be seen in Figure 1 
% in the project report.

% Forward difference formula for current1.dat file.
% The operations that I wrote in comments above are applied on the codes
% below.
for i = 1:8
    derivat_cur1_forward(i, 1) = (current1(i+1, 2) - current1(i, 2)) / delta_t1;
    voltage_e_cur1_forward(i, 1) = (derivat_cur1_forward(i) * L) + (R * current1(i, 2));
end

% Backward difference formula for current1.dat file.
for i = 2:8
    derivat_cur1_backward(i, 1) = (current1(i, 2) - current1(i-1, 2)) / delta_t1;
    voltage_e_cur1_backward(i, 1) = (derivat_cur1_backward(i) * L) + (R * current1(i, 2));
end

% Three point midpoint formula for current1.dat file.
for i = 2:8
    derivat_cur1_tpmp(i, 1)=((current1(i+1, 2) - current1(i-1, 2)) / (2*delta_t1));
    voltage_e_cur1_tpmp(i, 1)=(derivat_cur1_tpmp(i) * L) + (R*current1(i,2));
end

% Forward difference formula for current2.dat file.
for i = 1:12
    derivat_cur2_forward(i, 1) = (current2(i+1, 2) - current2(i, 2)) / delta_t2;
    voltage_e_cur2_forward(i, 1) = (derivat_cur2_forward(i) * L) + (R * current2(i, 2));
end

% Backward difference formula for current2.dat file.
for i = 2:12
    derivat_cur2_backward(i, 1) = (current2(i, 2) - current2(i-1, 2)) / delta_t2;
    voltage_e_cur2_backward(i, 1) = (derivat_cur2_backward(i) * L) + (R * current2(i, 2));
end

% Three point midpoint formula for current2.dat file.
for i = 2:12
    derivat_cur2_tpmp(i, 1)=(current2(i+1, 2) - current2(i-1, 2)) / (2*delta_t2);
    voltage_e_cur2_tpmp(i, 1)=(derivat_cur2_tpmp(i) * L) + (R*current2(i,2));
end

% Forward difference formula for current3.dat file.
for i = 1:24
    derivat_cur3_forward(i, 1) = (current3(i+1, 2) - current3(i, 2)) / delta_t3;
    voltage_e_cur3_forward(i, 1) = (derivat_cur3_forward(i) * L) + (R * current3(i, 2));
end

% Backward difference formula for current3.dat file.
for i = 2:24
    derivat_cur3_backward(i, 1) = (current3(i, 2) - current3(i-1, 2)) / delta_t3;
    voltage_e_cur3_backward(i, 1) = (derivat_cur3_backward(i) * L) + (R * current3(i, 2));
end

% Three point midpoint formula for current3.dat file.
for i = 2:24
    derivat_cur3_tpmp(i, 1)=(current3(i+1, 2) - current3(i-1, 2)) / (2*delta_t3);
    voltage_e_cur3_tpmp(i, 1)=(derivat_cur3_tpmp(i) * L) + (R * current3(i,2));
end

% Forward difference formula for current4.dat file.
for i = 1:60
    derivat_cur4_forward(i, 1) = (current4(i+1, 2) - current4(i, 2)) / delta_t4;
    voltage_e_cur4_forward(i, 1) = (derivat_cur4_forward(i) * L) + (R * current4(i, 2));
end

% Backward difference formula for current4.dat file.
for i = 2:60
    derivat_cur4_backward(i, 1) = (current4(i, 2) - current4(i-1, 2)) / delta_t4;
    voltage_e_cur4_backward(i, 1) = (derivat_cur4_backward(i) * L) + (R * current4(i, 2));
end

% Three point midpoint formula for current4.dat file.
for i = 2:60
    derivat_cur4_tpmp(i, 1)=(current4(i+1, 2) - current4(i-1, 2)) / (2*delta_t4);
    voltage_e_cur4_tpmp(i, 1)=(derivat_cur4_tpmp(i) * L) + (R * current4(i,2));
end

% In the datafiles, we have 9 rows in time section. To plot the data, I
% needed to delete last row of the time column which means I deleted the
% 0.6 ms and my plots does not show any value at the time 0.6 ms.
current1(end,:) = [];
current2(end,:) = [];
current3(end,:) = [];
current4(end,:) = [];

% I created new figures and plotted subplots with respect to appropriate
% data file and decorate it to make more understandable in the section
% below.
figure(1);
subplot(2,2,1)
plot(current1(:, 1), voltage_e_cur1_forward, '-*', current1(:, 1), voltage_e_cur1_backward, '-*', current1(:, 1), voltage_e_cur1_tpmp, '-*');
title('Voltage values of current1.dat');
xlabel('Time');
ylabel('Voltage');
legend('FD', 'BD', 'TPMP', 'Location', 'southeast');
legend('boxoff');

subplot(2,2,2)
plot(current2(:, 1), voltage_e_cur2_forward, '-*', current2(:, 1), voltage_e_cur2_backward, '-*', current2(:, 1), voltage_e_cur2_tpmp, '-*');
title('Voltage values of current2.dat');
xlabel('Time');
ylabel('Voltage');
legend('FD', 'BD', 'TPMP','Location','southeast');
legend('boxoff');

subplot(2,2,3)
plot(current3(:, 1), voltage_e_cur3_forward, '-*', current3(:, 1), voltage_e_cur3_backward, '-*', current3(:, 1), voltage_e_cur3_tpmp, '-*');
title('Voltage values of current3.dat');
xlabel('Time');
ylabel('Voltage');
legend('FD', 'BD', 'TPMP','Location','southeast');
legend('boxoff');

subplot(2,2,4)
plot(current4(:, 1), voltage_e_cur4_forward, '-*', current4(:, 1), voltage_e_cur4_backward, '-*', current4(:, 1), voltage_e_cur4_tpmp, '-*');
title('Voltage values of current4.dat');
xlabel('Time');
ylabel('Voltage');
legend('FD', 'BD', 'TPMP','Location','southeast');
legend('boxoff');

figure(2);
subplot(2,2,1)
plot(current1(:, 1), derivat_cur1_forward, '-*', current1(:, 1), derivat_cur1_backward, '-*', current1(:, 1), derivat_cur1_tpmp, '-*');
title('Derivative values of current1.dat');
xlabel('Time');
ylabel('Derivative');
legend('FD', 'BD', 'TPMP', 'Location', 'northeast');
legend('boxoff');

subplot(2,2,2)
plot(current2(:, 1), derivat_cur2_forward, '-*', current2(:, 1), derivat_cur2_backward, '-*', current2(:, 1), derivat_cur2_tpmp, '-*');
title('Derivative values of current2.dat');
xlabel('Time');
ylabel('Derivative');
legend('FD', 'BD', 'TPMP', 'Location', 'northeast');
legend('boxoff');

subplot(2,2,3)
plot(current3(:, 1), derivat_cur3_forward, '-*', current3(:, 1), derivat_cur3_backward, '-*', current3(:, 1), derivat_cur3_tpmp, '-*');
title('Derivative values of current3.dat');
xlabel('Time');
ylabel('Derivative');
legend('FD', 'BD', 'TPMP', 'Location', 'northeast');
legend('boxoff');

subplot(2,2,4)
plot(current4(:, 1), derivat_cur4_forward, '-*', current4(:, 1), derivat_cur4_backward, '-*', current4(:, 1), derivat_cur4_tpmp, '-*');
title('Derivative values of current4.dat');
xlabel('Time');
ylabel('Derivative');
legend('FD', 'BD', 'TPMP', 'Location', 'northeast');
legend('boxoff');

% Plot of currents with respect to time in the datafiles.
figure(3);
plot(current1(:,1), current1(:,2), current2(:,1), current2(:,2), current3(:,1), current3(:,2), current4(:,1), current4(:,2))
grid on;
title('Current values of current.dat files');
xlabel('Time');
ylabel('Currents');
legend('current1.dat', 'current2.dat', 'current3.dat', 'current4.dat', 'Location', 'southeast');
legend('boxoff');

% ERROR ANALYSIS
disp('***************ERRORS********************');
disp('The error for the current1.dat file;');
fprintf('When applied forward difference method = %f\n', (norm(voltage_e_cur1_forward - 12) / norm(voltage_e_cur1_forward)));
fprintf('When applied backward difference method = %f\n', (norm(voltage_e_cur1_backward - 12) / norm(voltage_e_cur1_backward)));
fprintf('When applied three point midpoint formula = %f\n', (norm(voltage_e_cur1_tpmp - 12) / norm(voltage_e_cur1_tpmp)));

fprintf('\nThe error for the current2.dat file;\n');
fprintf('When applied forward difference method = %f\n', (norm(voltage_e_cur2_forward - 12) / norm(voltage_e_cur2_forward)));
fprintf('When applied backward difference method = %f\n', (norm(voltage_e_cur2_backward - 12) / norm(voltage_e_cur2_backward)));
fprintf('When applied three point midpoint formula = %f\n', (norm(voltage_e_cur2_tpmp - 12) / norm(voltage_e_cur2_tpmp)));

fprintf('\nThe error for the current3.dat file;\n');
fprintf('When applied forward difference method = %f\n', (norm(voltage_e_cur3_forward - 12) / norm(voltage_e_cur3_forward)));
fprintf('When applied backward difference method = %f\n', (norm(voltage_e_cur3_backward - 12) / norm(voltage_e_cur3_backward)));
fprintf('When applied three point midpoint formula = %f\n', (norm(voltage_e_cur3_tpmp - 12) / norm(voltage_e_cur3_tpmp)));

fprintf('\nThe error for the current4.dat file;\n');
fprintf('When applied forward difference method = %f\n', (norm(voltage_e_cur4_forward - 12) / norm(voltage_e_cur4_forward)));
fprintf('When applied backward difference method = %f\n', (norm(voltage_e_cur4_backward - 12) / norm(voltage_e_cur4_backward)));
fprintf('When applied three point midpoint formula = %f\n', (norm(voltage_e_cur4_tpmp - 12) / norm(voltage_e_cur4_tpmp)));