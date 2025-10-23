% Validation script for IEEE 9-bus system
clear; clc; close all;

% Run the provided data file to load system parameters
ieee9_A2;

% Build admittance matrix from line data
nbus = 9;
Y = zeros(nbus, nbus);

% Add line admittances
for k = 1:length(nfrom)
    i = nfrom(k);
    j = nto(k);
    
    z = r(k) + 1j*x(k);
    y = 1/z;
    
    % Add series admittance
    Y(i,i) = Y(i,i) + y + 1j*b(k)/2;
    Y(j,j) = Y(j,j) + y + 1j*b(k)/2;
    Y(i,j) = Y(i,j) - y;
    Y(j,i) = Y(j,i) - y;
end

fprintf('=== IEEE 9-Bus System Power Flow Analysis ===\n');
fprintf('Admittance matrix built successfully.\n');

% Run Newton-Raphson power flow
[V, delta, Ps1, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display additional results
fprintf('\n=== Line Power Flows ===\n');
fprintf('From To    P(MW)    Q(Mvar)\n');
fprintf('---------------------------\n');

% Calculate line flows
for k = 1:length(nfrom)
    i = nfrom(k);
    j = nto(k);
    
    Vi = V(i) * exp(1j*delta(i));
    Vj = V(j) * exp(1j*delta(j));
    
    Iij = (Vi - Vj) * (1/(r(k) + 1j*x(k))) + Vi * (1j*b(k)/2);
    Sij = Vi * conj(Iij);
    
    fprintf('%2d  %2d   %7.2f   %7.2f\n', i, j, real(Sij)*Sbase, imag(Sij)*Sbase);
end

% Verify power balance
total_gen_P = Ps1 + sum(Pg(ipv)) + sum(Pg(ipq));
total_load_P = sum(Pd);
total_losses_P = total_gen_P - total_load_P;

fprintf('\n=== Power Balance ===\n');
fprintf('Total Generation: %.2f MW\n', total_gen_P);
fprintf('Total Load: %.2f MW\n', total_load_P);
fprintf('Total Losses: %.2f MW\n', total_losses_P);
