%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = grow(par,sol)            
            %% Set up.
            
            kgrid = par.kgrid; % Capital today (state variable).
            Agrid = par.Agrid; % Productivity (state variable).

            yout = sol.y; % Production function.
            kpol = sol.k; % Policy function for capital.
            cpol = sol.c; % Policy function for consumption.
            ipol = sol.i; % Policy function for investment.
            gpol = sol.g; % Policy function for government

            T = par.T; % Time periods.
            Asim = zeros(par.T*2,1); % Container for simulated productivity.
            ysim = zeros(par.T*2,1); % Container for simulated output.
            ksim = zeros(par.T*2,1); % Container for simulated capital stock.
            csim = zeros(par.T*2,1); % Container for simulated consumption.
            gsim = zeros(par.T*2,1); % Container for simulated government.
            isim = zeros(par.T*2,1); % Container for simulated investment.
            usim = zeros(par.T*2,1); % Container for simulated utility.

            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:); % Stationary distribution.
            cmat = cumsum(par.pmat,2); % CDF matrix.

            k0_ind = randsample(par.klen,1); % Index for initial capital stock.
            A0_ind = randsample(par.Alen,1,true,pmat0); % Index for initial productivity.

            Asim(1) = Agrid(A0_ind); % Productivity in period 1.
            ysim(1) = yout(k0_ind,A0_ind); % Output in period 1 given k0 and A0.
            csim(1) = cpol(k0_ind,A0_ind); % Consumption in period 1 given k0 and A0.
            gsim(1) = gpol(k0_ind,A0_ind); % Government in period 1 given k0 and A0.
            ksim(1) = kpol(k0_ind,A0_ind); % Capital choice for period 2 given k0 and A0.
            isim(1) = ipol(k0_ind,A0_ind); % Investment in period 1 given k0 and A0.
            usim(1) = model.utility(csim(1),gsim(1),par); % Utility in period 1 given k0 and A0.

            A1_ind = find(rand<=cmat(A0_ind,:)); % Draw productivity for next period.
            A0_ind = A1_ind(1);

            %% Simulate endogenous and exogenous variables.

            for j = 2:T*2 % Time loop.
                kt_ind = find(ksim(j-1)==kgrid); % Capital choice in the previous period is the state today. Find where the latter is on the grid.
                Asim(j) = Agrid(A0_ind); % Productivity in period t.
                ysim(j) = yout(kt_ind,A0_ind); % Output in period t.
                csim(j) = cpol(kt_ind,A0_ind); % Consumption in period t.
                gsim(j) = gpol(kt_ind,A0_ind); % Government in period t.
                ksim(j) = kpol(kt_ind,A0_ind); % Capital stock for period t+1.
                isim(j) = ipol(kt_ind,A0_ind); % Investment in period t.
                usim(j) = model.utility(csim(j),gsim(j),par); % Utility in period t.
                A1_ind = find(rand<=cmat(A0_ind,:)); % Draw next state.
                A0_ind = A1_ind(1); % State next period.
            end

            sim = struct();
            
            % Burn the first half.
            sim.Asim = Asim(T+1:2*T,1); % Simulated productivity.
            sim.ysim = ysim(T+1:2*T,1); % Simulated output.
            sim.ksim = ksim(T+1:2*T,1); % Simulated capital choice.
            sim.csim = csim(T+1:2*T,1); % Simulated consumption.
            sim.gsim = gsim(T+1:2*T,1); % Simulated government.
            sim.isim = isim(T+1:2*T,1); % Simulated investment.
            sim.usim = usim(T+1:2*T,1); % Simulated utility.
             
        end
        
    end
end