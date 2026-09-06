% This script explores some of the parameter space of the anticipatory
% frequency adaptive delayed self-feedback Kuramoto-style model introduced 
% in the manuscript "Delay and amplification of auditory feedback for 
% walking: effects on variability, cadence, and temporal cortex activity".
% Details and more help can be found in anticipatory_kuramoto_rk4.m. 
% It uses an RK4 solver to run the model with the given initial conditions.
% Running the master script will explore a grid on the parameter space 
% of k_2 and tau_2, and it will repeat the simulation the given number of
% times for each combination of parameters. Because of iterations, running
% the entire loop in the master script will take a few hours on a 
% consumer-grade computer.
% Dobri Dotov, UNO, 2026

REZ = table('Size', [0, 11], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'k1','k2','epsilon1','epsilon2','taum','taus1','taus2','relphase0','sigma','deltaOmega','CV'});

repetitions = 1e2; % For the paper this was set to 5e2.
for n =  1:repetitions
    % Parameter space
    k_follower_1 = 2;
    k_follower_2_vec = linspace(0,4,21);
    epsilon_1 = k_follower_1*2;
    tau_leader = pi/6;
    tau_follower_1 = pi/5;
    tau_follower_2_vec = 0:pi/10:2*pi;
    init_rel_phase = 0;
    sigma = 1;
    [X, Y] = ndgrid(tau_follower_2_vec,k_follower_2_vec);
    param_space = [X(:) Y(:)];
    tau_follower_2_vec = param_space(:,1);
    k_follower_2_vec = param_space(:,2);

    % Initialize variables
    Z = zeros(size(param_space,1),1);
    cycle_variability = zeros(size(param_space,1),1);
    rel_phase_r = zeros(size(param_space,1),1);
    rel_phase = zeros(size(param_space,1),1);
    rel_phase_error = zeros(size(param_space,1),1);
    tau = zeros(size(param_space,1),2);
    omega_delta = zeros(size(param_space,1),1);

    for r = 1:size(param_space,1)
        K = [k_follower_1 k_follower_2_vec(r,1)];
        epsilon_2 =  k_follower_2_vec(r,1)*2;
        epsilon = [epsilon_1 epsilon_2];
        tau_arg = [tau_leader tau_follower_1 tau_follower_2_vec(r,1)];
        [cycle_variability(r),rel_phase_r(r),rel_phase(r),rel_phase_error(r),tau(r,:),omega_delta(r)] = ...
            anticipatory_kuramoto_rk4(K,epsilon,tau_arg,init_rel_phase,sigma,0,0);
    end

    rez = array2table([k_follower_1.*ones(size(k_follower_2_vec)) k_follower_2_vec epsilon_1.*ones(size(k_follower_2_vec)) k_follower_2_vec*2 ...
        tau_leader.*ones(size(k_follower_2_vec)) tau_follower_1.*ones(size(k_follower_2_vec)) tau_follower_2_vec ...
        init_rel_phase.*ones(size(k_follower_2_vec)) sigma.*ones(size(k_follower_2_vec)) ...
        omega_delta cycle_variability*1e2], 'VariableNames', REZ.Properties.VariableNames);
    REZ = [REZ; rez];

    fprintf('%6.2f%%', n/repetitions*1e2);
    if mod(n,10)==0; fprintf('\n'); end
end
% writetable(REZ,['rez_' char(datetime("now",'Format','yyyy-MM-dd')) '.csv']);

% Visualize
figures_across_runs;