% validate_ieee9.m - Validation script for IEEE 9-bus system
% This script validates the Newton-Raphson power flow implementation
% using the IEEE 9-bus test system.

clear; clc; close all;

fprintf('=== IEEE 9-Bus System Power Flow Validation ===\n\n');

% Load the IEEE 9-bus data
fprintf('Loading IEEE 9-bus system data...\n');
ieee9_A2;

% Build admittance matrix
fprintf('Building admittance matrix...\n');
nbus = 9;
Y = zeros(nbus, nbus);

for k = 1:length(nfrom)
    i = nfrom(k);
    j = nto(k);
    
    z = r(k) + 1j*x(k);
    y_series = 1/z;
    y_shunt = 1j*b(k)/2;
    
    Y(i,j) = Y(i,j) - y_series;
    Y(j,i) = Y(j,i) - y_series;
    Y(i,i) = Y(i,i) + y_series + y_shunt;
    Y(j,j) = Y(j,j) + y_series + y_shunt;
end

% Display system information
fprintf('System Information:\n');
fprintf('  Number of buses: %d\n', nbus);
fprintf('  Slack bus: %d\n', is);
fprintf('  PV buses: ['); fprintf('%d ', ipv); fprintf(']\n');
fprintf('  PQ buses: ['); fprintf('%d ', ipq); fprintf(']\n');
fprintf('  Power base: %.0f MVA\n', Sbase);
fprintf('  Tolerance: %.0e p.u.\n', toler);
fprintf('  Max iterations: %d\n\n', maxiter);

% Run power flow
fprintf('Running Newton-Raphson power flow...\n');
[V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display results
fprintf('\n=== POWER FLOW RESULTS ===\n');
fprintf('Convergence: %d iterations, %.4f seconds\n\n', N, time);

fprintf('Bus Results:\n');
fprintf('Bus  Type    V(p.u.)   Angle(deg)   Pg(MW)   Qg(Mvar)   Pd(MW)   Qd(Mvar)\n');
fprintf('---- ------ --------- ----------- --------- ---------- --------- ----------\n');

for i = 1:nbus
    angle_deg = rad2deg(delta(i));
    
    if i == is
        bus_type = 'Slack';
    elseif ismember(i, ipv)
        bus_type = 'PV';
    else
        bus_type = 'PQ';
    end
    
    fprintf('%3d  %-6s %9.4f %11.4f %9.1f %10.1f %9.1f %9.1f\n', ...
            i, bus_type, V(i), angle_deg, Pg(i), Qg(i), Pd(i), Qd(i));
end

fprintf('\nPV Bus Reactive Generation:\n');
for i = 1:length(ipv)
    fprintf('  Bus %d: Qg = %.2f Mvar\n', ipv(i), Qgv(i));
end

fprintf('\nSlack Bus Generation: P = %.2f MW\n', Psl);

% Power balance check
P_total_gen = sum(Pg) + Psl;
P_total_load = sum(Pd);
P_loss = P_total_gen - P_total_load;

fprintf('\nPower Balance:\n');
fprintf('  Total Generation: %.2f MW\n', P_total_gen);
fprintf('  Total Load:      %.2f MW\n', P_total_load);
fprintf('  Total Losses:    %.2f MW\n', P_loss);

% Plot results
figure('Position', [100, 100, 1000, 700]);

subplot(2,2,1);
bar(V);
title('Voltage Magnitude Profile');
xlabel('Bus Number');
ylabel('Voltage (p.u.)');
grid on;
ylim([0.95, 1.05]);

subplot(2,2,2);
bar(rad2deg(delta));
title('Voltage Angle Profile');
xlabel('Bus Number');
ylabel('Angle (degrees)');
grid on;

subplot(2,2,3);
% Bus type distribution
bus_types = {'Slack', 'PV', 'PQ'};
type_counts = [1, length(ipv), length(ipq)];
pie(type_counts, bus_types);
title('Bus Type Distribution');

subplot(2,2,4);
% Power summary
power_data = [P_total_gen, P_total_load, P_loss; ...
              sum(Qg) + sum(Qgv), sum(Qd), sum(Qg) + sum(Qgv) - sum(Qd)];
bar(power_data);
title('Power Summary');
set(gca, 'XTickLabel', {'Active Power', 'Reactive Power'});
ylabel('Power (MW/Mvar)');
legend('Generation', 'Load', 'Losses', 'Location', 'northeast');
grid on;

fprintf('\nValidation completed successfully!\n');
