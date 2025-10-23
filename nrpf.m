function [V, delta, Ps1, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%NRPF Newton-Raphson Power Flow Solver
%   [V, delta, Ps1, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%
%   Inputs:
%   Y - Network admittance matrix
%   is - Index of slack node
%   ipq - Vector of PQ node indices
%   ipv - Vector of PV node indices  
%   Pg, Qg - Generation active/reactive power vectors (MW/Mvar)
%   Pd, Qd - Demand active/reactive power vectors (MW/Mvar)
%   V0 - Initial voltage magnitude vector (p.u.)
%   Sbase - Network power base (MVA)
%   toler - Convergence tolerance (p.u.)
%   maxiter - Maximum iterations
%
%   Outputs:
%   V - Voltage magnitudes (p.u.)
%   delta - Voltage angles (radians)
%   Ps1 - Slack bus active generation (MW)
%   Qgv - PV bus reactive generation (Mvar)
%   N - Number of iterations to convergence
%   time - CPU time (seconds)

    tstart = tic; % Start timing
    
    % Convert power inputs to per unit
    Pg_pu = Pg / Sbase;
    Qg_pu = Qg / Sbase;
    Pd_pu = Pd / Sbase;
    Qd_pu = Qd / Sbase;
    
    nbus = size(Y, 1); % Total number of buses
    
    % Initialize voltage vectors
    V = V0; % Voltage magnitudes
    delta = zeros(nbus, 1); % Voltage angles (start flat)
    
    % Combine all non-slack buses
    non_slack = sort([ipq; ipv]);
    
    % Separate PQ and PV buses for mismatch calculations
    npq = length(ipq);
    npv = length(ipv);
    
    fprintf('Starting Newton-Raphson Power Flow...\n');
    fprintf('System: %d buses (%d PQ, %d PV, 1 Slack)\n', nbus, npq, npv);
    fprintf('Tolerance: %.2e, Max iterations: %d\n', toler, maxiter);
    
    for N = 1:maxiter
        % Calculate power mismatches
        [dP, dQ, Pcalc, Qcalc] = calculate_mismatches(Y, V, delta, Pg_pu, Qg_pu, Pd_pu, Qd_pu, ipq, ipv, is);
        
        % Check convergence
        max_mismatch = max([abs(dP); abs(dQ)]);
        fprintf('Iteration %d: Max mismatch = %.6f p.u.\n', N, max_mismatch);
        
        if max_mismatch < toler
            fprintf('Converged in %d iterations!\n', N);
            break;
        end
        
        % Build Jacobian matrix
        J = build_jacobian(Y, V, delta, ipq, ipv, non_slack);
        
        % Build mismatch vector for non-slack buses
        mismatch = [dP(non_slack); dQ(ipq)];
        
        % Solve for corrections using linsolve
        corrections = linsolve(J, -mismatch);
        
        % Extract angle and voltage corrections
        n_non_slack = length(non_slack);
        ddelta_non_slack = corrections(1:n_non_slack);
        dV_pq = corrections(n_non_slack+1:end);
        
        % Update voltage angles for non-slack buses
        delta(non_slack) = delta(non_slack) + ddelta_non_slack;
        
        % Update voltage magnitudes for PQ buses only
        V(ipq) = V(ipq) + dV_pq;
        
        % Ensure PV bus voltages remain at specified values
        V(ipv) = V0(ipv);
    end
    
    if N == maxiter && max_mismatch >= toler
        warning('Power flow did not converge within maximum iterations!');
    end
    
    % Calculate final power injections
    [~, ~, Pcalc, Qcalc] = calculate_mismatches(Y, V, delta, Pg_pu, Qg_pu, Pd_pu, Qd_pu, ipq, ipv, is);
    
    % Calculate output quantities
    Ps1 = Pcalc(is) * Sbase; % Slack bus active power in MW
    Qgv = Qcalc(ipv) * Sbase; % PV bus reactive power in Mvar
    
    time = toc(tstart);
    
    fprintf('\n=== Power Flow Results ===\n');
    fprintf('Solution time: %.4f seconds\n', time);
    fprintf('Iterations: %d\n', N);
    fprintf('Slack bus power: %.3f MW\n', Ps1);
    
    % Display bus results
    fprintf('\nBus Results:\n');
    fprintf('Bus  Type  |V|(p.u.)  Angle(deg)   Pg(MW)   Qg(Mvar)   Pd(MW)   Qd(Mvar)\n');
    fprintf('------------------------------------------------------------------------\n');
    
    for i = 1:nbus
        angle_deg = rad2deg(delta(i));
        if i == is
            bus_type = 'Slack';
            Pg_disp = Ps1;
            Qg_disp = Qcalc(i) * Sbase;
        elseif ismember(i, ipv)
            bus_type = 'PV';
            Pg_disp = Pg(i);
            Qg_disp = Qgv(ipv == i);
        else
            bus_type = 'PQ';
            Pg_disp = Pg(i);
            Qg_disp = Qg(i);
        end
        fprintf('%2d   %-5s  %8.4f   %8.2f   %8.1f   %8.1f   %8.1f   %8.1f\n', ...
            i, bus_type, V(i), angle_deg, Pg_disp, Qg_disp, Pd(i), Qd(i));
    end
end

function [dP, dQ, Pcalc, Qcalc] = calculate_mismatches(Y, V, delta, Pg_pu, Qg_pu, Pd_pu, Qd_pu, ipq, ipv, is)
% Calculate power mismatches and computed power injections
    nbus = length(V);
    Pcalc = zeros(nbus, 1);
    Qcalc = zeros(nbus, 1);
    
    % Calculate computed power injections
    for i = 1:nbus
        for k = 1:nbus
            theta_ik = delta(i) - delta(k);
            Pcalc(i) = Pcalc(i) + V(i) * V(k) * (real(Y(i,k)) * cos(theta_ik) + imag(Y(i,k)) * sin(theta_ik));
            Qcalc(i) = Qcalc(i) + V(i) * V(k) * (real(Y(i,k)) * sin(theta_ik) - imag(Y(i,k)) * cos(theta_ik));
        end
    end
    
    % Calculate mismatches
    dP = Pg_pu - Pd_pu - Pcalc;
    dQ = Qg_pu - Qd_pu - Qcalc;
    
    % Zero out slack bus mismatches (they're not in the Jacobian)
    dP(is) = 0;
    
    % Zero out PV bus reactive power mismatches
    dQ(ipv) = 0;
end

function J = build_jacobian(Y, V, delta, ipq, ipv, non_slack)
% Build the Jacobian matrix
    nbus = length(V);
    npq = length(ipq);
    n_non_slack = length(non_slack);
    
    % Initialize Jacobian
    J = zeros(n_non_slack + npq, n_non_slack + npq);
    
    % Fill J11 (dP/ddelta)
    for i = 1:n_non_slack
        row_idx = i;
        bus_i = non_slack(i);
        for j = 1:n_non_slack
            col_idx = j;
            bus_j = non_slack(j);
            
            if bus_i == bus_j
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    if k ~= bus_i
                        theta_ik = delta(bus_i) - delta(k);
                        sum_term = sum_term + V(bus_i) * V(k) * ...
                            (-real(Y(bus_i,k)) * sin(theta_ik) + imag(Y(bus_i,k)) * cos(theta_ik));
                    end
                end
                J(row_idx, col_idx) = sum_term;
            else
                % Off-diagonal element
                theta_ij = delta(bus_i) - delta(bus_j);
                J(row_idx, col_idx) = V(bus_i) * V(bus_j) * ...
                    (real(Y(bus_i,bus_j)) * sin(theta_ij) - imag(Y(bus_i,bus_j)) * cos(theta_ij));
            end
        end
    end
    
    % Fill J12 (dP/dV)
    for i = 1:n_non_slack
        row_idx = i;
        bus_i = non_slack(i);
        for j = 1:npq
            col_idx = n_non_slack + j;
            bus_j = ipq(j);
            
            if bus_i == bus_j
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    theta_ik = delta(bus_i) - delta(k);
                    sum_term = sum_term + V(k) * (real(Y(bus_i,k)) * cos(theta_ik) + imag(Y(bus_i,k)) * sin(theta_ik));
                end
                J(row_idx, col_idx) = V(bus_i) * sum_term + V(bus_i) * real(Y(bus_i,bus_i)) * cos(0);
            else
                % Off-diagonal element
                theta_ij = delta(bus_i) - delta(bus_j);
                J(row_idx, col_idx) = V(bus_i) * (real(Y(bus_i,bus_j)) * cos(theta_ij) + imag(Y(bus_i,bus_j)) * sin(theta_ij));
            end
        end
    end
    
    % Fill J21 (dQ/ddelta)
    for i = 1:npq
        row_idx = n_non_slack + i;
        bus_i = ipq(i);
        for j = 1:n_non_slack
            col_idx = j;
            bus_j = non_slack(j);
            
            if bus_i == bus_j
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    if k ~= bus_i
                        theta_ik = delta(bus_i) - delta(k);
                        sum_term = sum_term + V(bus_i) * V(k) * ...
                            (real(Y(bus_i,k)) * cos(theta_ik) + imag(Y(bus_i,k)) * sin(theta_ik));
                    end
                end
                J(row_idx, col_idx) = sum_term;
            else
                % Off-diagonal element
                theta_ij = delta(bus_i) - delta(bus_j);
                J(row_idx, col_idx) = -V(bus_i) * V(bus_j) * ...
                    (real(Y(bus_i,bus_j)) * cos(theta_ij) + imag(Y(bus_i,bus_j)) * sin(theta_ij));
            end
        end
    end
    
    % Fill J22 (dQ/dV)
    for i = 1:npq
        row_idx = n_non_slack + i;
        bus_i = ipq(i);
        for j = 1:npq
            col_idx = n_non_slack + j;
            bus_j = ipq(j);
            
            if bus_i == bus_j
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    theta_ik = delta(bus_i) - delta(k);
                    sum_term = sum_term + V(k) * (real(Y(bus_i,k)) * sin(theta_ik) - imag(Y(bus_i,k)) * cos(theta_ik));
                end
                J(row_idx, col_idx) = V(bus_i) * sum_term - V(bus_i) * imag(Y(bus_i,bus_i)) * cos(0);
            else
                % Off-diagonal element
                theta_ij = delta(bus_i) - delta(bus_j);
                J(row_idx, col_idx) = V(bus_i) * (real(Y(bus_i,bus_j)) * sin(theta_ij) - imag(Y(bus_i,bus_j)) * cos(theta_ij));
            end
        end
    end
end
