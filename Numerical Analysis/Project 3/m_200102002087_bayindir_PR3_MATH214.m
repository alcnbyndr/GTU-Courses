close all; clear all; clc;

L = 0.1; % 100 mH = 0.1 H
DT = 0.025;
stored_energy_simpsons(1,41) = zeros;
stored_energy_midpoint(1,41) = zeros;
stored_energy_trapezoidal(1,41) = zeros;

load pr3data.dat;

time_values = pr3data(:, 1)';
current_values = pr3data(:, 2)';
voltage_values = pr3data(:, 3)';

for i = 2:length(time_values)
    derivat_cur1_backward(i) = (current_values(i) - current_values(i-1)) / DT;
    voltage_e_cur1_backward(i) = (derivat_cur1_backward(i) * L);
end

% COmposite Simpson's rule
for j = 1:length(time_values) 
    power(j) = voltage_values(j) * current_values(j);
    for k = 1:j
        if (k==1 || k==i)
            stored_energy_simpsons(j) = stored_energy_simpsons(j) + power(k);
        elseif(rem(k, 2) == 0)
            stored_energy_simpsons(j) = stored_energy_simpsons(j) + 4*power(k);
        elseif(rem(k, 2) == 1)
            stored_energy_simpsons(j) = stored_energy_simpsons(j) + 2*power(k);
        end
    end
    stored_energy_simpsons(j) = stored_energy_simpsons(j) * DT / 3;
end

% Composite Midpoint
for x = 1:length(time_values) - 1
    for t = 1:2:x
        stored_energy_midpoint(x) = stored_energy_midpoint(x) + 2*DT*power(t);
    end
end

for p = 2:length(time_values) - 1
    n = time_values(p)/DT;
    h = (time_values(p) - time_values(1)) / n; 
    for y = 1:(n+1)
        if (j==1 || j==(n+1)) 
            stored_energy_trapezoidal(p) = stored_energy_trapezoidal(p) + power(y);
        else
            stored_energy_trapezoidal(p) = stored_energy_trapezoidal(p) + 2 * power(y);
        end
    end
    stored_energy_trapezoidal(p) = (h/2) * stored_energy_trapezoidal(p);
end

voltage_values_conc = vertcat(time_values, voltage_values, voltage_e_cur1_backward)';
 
% For the fourth question
stored_energy_eq3 = (L .* (current_values).^2) ./ 2;

plot(time_values(1,2:end), current_values(1,2:end), time_values(1,2:end), ...
voltage_values(1,2:end),'k', time_values(1,2:end), voltage_e_cur1_backward(1,2:end), 'r--*', 'LineWidth', 2);
grid on;
title('Current/Voltage - Time Graph');
xlabel('Time'); ylabel('Current (I) AND Voltage (V)');
legend('Current', 'Voltage', 'Converged Voltage');

figure(2);
plot(time_values, stored_energy_eq3, '*', time_values, stored_energy_simpsons, time_values, stored_energy_midpoint, time_values, stored_energy_trapezoidal, 'LineWidth', 2);
grid on;
title('Results of methods');
xlabel('Time'); ylabel('Values');
legend('Stored Energy in real', 'Composite Simpson Method', 'Composite Midpoint Method', 'Composite Trapezoidal Method');

% Error Calculations
fprintf('Error rate of composite Simpson formula = %f\n', (norm(stored_energy_simpsons - stored_energy_eq3) / norm(stored_energy_simpsons)));
fprintf('Error rate of composite midpoint formula = %f\n', (norm(stored_energy_midpoint - stored_energy_eq3) / norm(stored_energy_midpoint)));
fprintf('Error rate of composite trapezoidal method = %f\n', (norm(stored_energy_trapezoidal - stored_energy_eq3) / norm(stored_energy_trapezoidal)));
