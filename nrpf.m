function [V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%NRPF Newton-Raphson Power Flow Solution
%   [V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%   solves the power flow problem using the Newton-Raphson method.

    tstart = tic; % Start timer
    
    % Convert inputs to per unit
    Pg_pu = Pg / Sbase;
    Qg_pu = Qg / Sbase;
    Pd_pu = Pd / Sbase;
    Qd_pu = Qd / Sbase;
    
    nbus = size(Y, 1); % Total number of buses
    
    % Initialize variables
    V = V0; % Voltage magnitudes
    delta = zeros(nbus, 1); % Voltage angles
    P_spec = Pg_pu - Pd_pu; % Net specified active power
    Q_spec = Qg_pu - Qd_pu; % Net specified reactive power
    
    % Identify all buses and create index mapping
    all_buses = 1:nbus;
    pq_buses = ipq;
    pv_buses = ipv;
    
    % Create mismatch vectors indices
    npq = length(pq_buses);
    npv = length(pv_buses);
    
    % Total number of unknowns
    nunknowns = 2*npq + npv;
    
    % Build index mapping for state variables
    delta_idx = zeros(nbus, 1);
    V_idx = zeros(nbus, 1);
    
    idx = 1;
    % PV buses first (only delta unknown)
    for i = 1:npv
        bus = pv_buses(i);
        delta_idx(bus) = idx;
        idx = idx + 1;
    end
    
    % PQ buses (both delta and V unknown)
    for i = 1:npq
        bus = pq_buses(i);
        delta_idx(bus) = idx;
        idx = idx + 1;
        V_idx(bus) = idx;
        idx = idx + 1;
    end
    
    % Newton-Raphson iteration
    converged = false;
    N = 0;
    
    while ~converged && N < maxiter
        N = N + 1;
        
        % Calculate power injections
        P_calc = zeros(nbus, 1);
        Q_calc = zeros(nbus, 1);
        
        for i = 1:nbus
            theta_i = delta(i);
            V_i = V(i);
            
            for j = 1:nbus
                theta_ij = theta_i - delta(j);
                Y_mag = abs(Y(i,j));
                Y_angle = angle(Y(i,j));
                
                P_calc(i) = P_calc(i) + V_i * V(j) * Y_mag * cos(theta_ij + Y_angle);
                Q_calc(i) = Q_calc(i) + V_i * V(j) * Y_mag * sin(theta_ij + Y_angle);
            end
        end
        
        % Calculate mismatches
        deltaP = P_spec - P_calc;
        deltaQ = Q_spec - Q_calc;
        
        % Remove slack bus from mismatch
        deltaP(is) = [];
        
        % Create Q mismatch only for PQ buses
        deltaQ_pq = [];
        for i = 1:nbus
            if i ~= is && ismember(i, pq_buses)
                deltaQ_pq = [deltaQ_pq; deltaQ(i)];
            end
        end
        
        % Create complete mismatch vector
        mismatch = [deltaP; deltaQ_pq];
        
        % Check convergence
        max_mismatch = max(abs(mismatch));
        if max_mismatch < toler
            converged = true;
            break;
        end
        
        % Build Jacobian matrix
        J = zeros(nunknowns, nunknowns);
        
        % 1. dP/ddelta terms
        for i = 1:nbus
            if i == is, continue; end
            
            row_idx = delta_idx(i);
            if row_idx == 0, continue; end
            
            for j = 1:nbus
                if j == is, continue; end
                
                theta_ij = delta(i) - delta(j);
                Y_mag = abs(Y(i,j));
                Y_angle = angle(Y(i,j));
                
                if i == j
                    % Diagonal elements
                    J(row_idx, row_idx) = -Q_calc(i) - (V(i)^2) * imag(Y(i,i));
                else
                    % Off-diagonal elements
                    col_idx = delta_idx(j);
                    if col_idx > 0
                        J(row_idx, col_idx) = V(i) * V(j) * Y_mag * sin(theta_ij + Y_angle);
                    end
                end
            end
        end
        
        % 2. dP/dV terms
        for i = 1:nbus
            if i == is, continue; end
            
            row_idx = delta_idx(i);
            if row_idx == 0, continue; end
            
            for j = 1:nbus
                if ~ismember(j, pq_buses), continue; end
                
                theta_ij = delta(i) - delta(j);
                Y_mag = abs(Y(i,j));
                Y_angle = angle(Y(i,j));
                
                col_idx = V_idx(j);
                if col_idx == 0, continue; end
                
                if i == j
                    J(row_idx, col_idx) = P_calc(i)/V(i) + V(i) * real(Y(i,i));
                else
                    J(row_idx, col_idx) = V(i) * Y_mag * cos(theta_ij + Y_angle);
                end
            end
        end
        
        % 3. dQ/ddelta terms
        for i = 1:nbus
            if i == is || ~ismember(i, pq_buses), continue; end
            
            row_idx = V_idx(i);
            if row_idx == 0, continue; end
            
            for j = 1:nbus
                if j == is, continue; end
                
                theta_ij = delta(i) - delta(j);
                Y_mag = abs(Y(i,j));
                Y_angle = angle(Y(i,j));
                
                if i == j
                    col_idx = delta_idx(i);
                    if col_idx > 0
                        J(row_idx, col_idx) = P_calc(i) - (V(i)^2) * real(Y(i,i));
                    end
                else
                    col_idx = delta_idx(j);
                    if col_idx > 0
                        J(row_idx, col_idx) = -V(i) * V(j) * Y_mag * cos(theta_ij + Y_angle);
                    end
                end
            end
        end
        
        % 4. dQ/dV terms
        for i = 1:nbus
            if i == is || ~ismember(i, pq_buses), continue; end
            
            row_idx = V_idx(i);
            if row_idx == 0, continue; end
            
            for j = 1:nbus
                if ~ismember(j, pq_buses), continue; end
                
                theta_ij = delta(i) - delta(j);
                Y_mag = abs(Y(i,j));
                Y_angle = angle(Y(i,j));
                
                col_idx = V_idx(j);
                if col_idx == 0, continue; end
                
                if i == j
                    J(row_idx, col_idx) = Q_calc(i)/V(i) - V(i) * imag(Y(i,i));
                else
                    J(row_idx, col_idx) = V(i) * Y_mag * sin(theta_ij + Y_angle);
                end
            end
        end
        
        % Solve using linsolve()
        opts.SYM = false;
        opts.POSDEF = false;
        correction = linsolve(J, mismatch, opts);
        
        % Update state variables
        for i = 1:nbus
            if i == is, continue; end
            idx_delta = delta_idx(i);
            if idx_delta > 0
                delta(i) = delta(i) + correction(idx_delta);
            end
        end
        
        for i = 1:npq
            bus = pq_buses(i);
            idx_V = V_idx(bus);
            if idx_V > 0
                V(bus) = V(bus) + correction(idx_V);
            end
        end
    end
    
    % Calculate final power injections
    P_calc_final = zeros(nbus, 1);
    Q_calc_final = zeros(nbus, 1);
    
    for i = 1:nbus
        for j = 1:nbus
            theta_ij = delta(i) - delta(j);
            Y_mag = abs(Y(i,j));
            Y_angle = angle(Y(i,j));
            
            P_calc_final(i) = P_calc_final(i) + V(i) * V(j) * Y_mag * cos(theta_ij + Y_angle);
            Q_calc_final(i) = Q_calc_final(i) + V(i) * V(j) * Y_mag * sin(theta_ij + Y_angle);
        end
    end
    
    % Output results
    Psl = P_calc_final(is) * Sbase;
    
    Qgv = zeros(length(ipv), 1);
    for i = 1:length(ipv)
        bus = ipv(i);
        Qgv(i) = Q_calc_final(bus) * Sbase;
    end
    
    time = toc(tstart);
    
    if ~converged
        warning('Newton-Raphson failed to converge in %d iterations. Max mismatch: %e', maxiter, max_mismatch);
    else
        fprintf('Newton-Raphson converged in %d iterations\n', N);
    end
end
