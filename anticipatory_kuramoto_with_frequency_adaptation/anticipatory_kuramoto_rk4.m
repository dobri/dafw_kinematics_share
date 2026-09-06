function [cycle_variability,mean_rel_phase_consistency,mean_rel_phase,rel_phase_error,tau_s,omega_delta,omega_s] = anticipatory_kuramoto_rk4(varargin)
%ANTICIPATORY_KURAMOTO_RK4  Anticipatory Kuramoto model with multiple delayed feedback streams and frequency adaptation
%
%       [cycle_variability,mean_rel_phase_consistency,mean_rel_phase,rel_phase_error,tau_s,omega_delta,omega_s] = ...
%                                       anticipatory_kuramoto_rk4(varargin)
%
%       Arguments:
%       [k_1 and k_2]
%       [epsilon_1 epsilon_2]
%       [tau_leader tau_follower_1 tau_follower_2]
%       relative_phase_initial of follower relative to leader's zero phase.
%       show figs flag
%       print figs flag
%
%Dobri Dotov, UNO, 2026

if numel(varargin) > 0
    k1       = varargin{1}(1);   % Phase coupling strength
    k2       = varargin{1}(2);
    eps1     = varargin{2}(1);   % Frequency adaptation strength
    eps2     = varargin{2}(2);
    tau_m    = varargin{3}(1);
    tau_s    = varargin{3}(2:3); % Delay phase (anticipation horizon)
    theta_m0 = 0;
    theta_s0 = -varargin{4}; % Initial rel phase
    sigma    = varargin{5};
    plotting_flag = varargin{6}; % Viz
    save_fig = varargin{7};
else % Default parameters
    k1       = 2;
    k2       = 2;
    eps1     = 5;
    eps2     = 5;
    tau_m    = 2*pi*.2;
    tau_s    = [2*2*pi 2*2*pi];
    sigma    = .0;
    theta_m0 = 0;
    theta_s0 = 0;
    plotting_flag = 1;
    save_fig = 0;
end

%% Natural frequencies
% Leader doesn't have a natural frequency. Its frequency is the follower's.
omega_s0 = 2*2*pi;                  % Follower freq 0
tau_L    = tau_m/omega_s0(1);       % leader delay (s)
tau_F1   = tau_s(1)/omega_s0(1);    % follower delay 1 (s)
tau_F2   = tau_s(2)/omega_s0(1);
mu       = 10;                      % freq adaptation timescale
dt       = 0.001;                   % Time step
T        = 120;                     % Total time
n        = round(T/dt+1);
alpha    = 0;                       % Ignore for now

%% Run
[t, theta_m, theta_s, omega_s, coupling_term] = run_sim(omega_s0, theta_m0, theta_s0, ...
    k1, k2, eps1, eps2, mu, alpha, tau_L, tau_F1, tau_F2, sigma, dt, n);
omega_m = omega_s;


%% Visualize
% Calculate synchronization measures
% Wrapped phase difference for plotting
% Phase difference
rel_phase_wrapped = wrap_angle(theta_m - theta_s);
z = exp(1i * rel_phase_wrapped);
rel_phase = angle(z);
mean_rel_phase = angle(mean(z(round(numel(z)/2):end)));
mean_rel_phase_consistency = abs(mean(z(round(numel(z)/2):end)));
final_rel_phase = angle(mean(z(end-round(1/dt):end)));
omega_delta = (omega_s(end) - omega_s(1))/omega_s(1)*1e2;
rel_phase_error = mod(mean_rel_phase - tau_m + tau_s(1)*k1/(k1+k2) + tau_s(2)*k2/(k1+k2) + pi, 2*pi) - pi;

% CV in steady state
phase_wrapped = mod(theta_s(round(numel(theta_s)/2):end), 2*pi);
increments = find(diff(phase_wrapped)<-pi)+1;
step_times = [];
previous_step = -inf;
for inc = increments
    if t(inc) > (previous_step+.2)
        step_times = vertcat(step_times, t(inc)); %#ok<AGROW>
        previous_step = t(inc);
    end
end
cycle_variability = std(diff(step_times))/mean(diff(step_times));


%% Plotting
if plotting_flag == 1
    % Display results
    fprintf('External delay τ_m: %.2f rad\n', tau_m);
    fprintf('Self-delay τ_1 and τ_2: %.2f and %.2f rad\n', tau_s);
    fprintf('Phase coupling K_1 and K_2: %.2f and %.2f\n', [k1 k2]);
    fprintf('Freq coupling ε_1 and ε_2: %.2f %.2f\n', [eps1 eps2]);
    fprintf('Final phase lag: %.6f rad\n', final_rel_phase);
    fprintf('Change in omega: %.6f%%\n', omega_delta);
    fprintf('Cycle variability: %.6f%%\n', cycle_variability*1e2);
    fprintf('\n');

    % Phase evolution
    figure(25421)
    subplot(3, 1, 1)
    plot(t, mod(theta_m,2*pi), 'k-', 'LineWidth', 2, 'DisplayName', 'Leader')
    hold on
    plot(t, mod(theta_s,2*pi), 'm-', 'LineWidth', 2, 'DisplayName', 'Follower')
    hold off
    xlabel('Time, s')
    ylabel('\theta [degrees]')
    legend('show','location','north')
    grid off
    set(gca,'YTick',[0 pi 2*pi],'YtickLabel',[0 pi 2*pi]./2/pi*360)
    ylim([0 2*pi])
    set(gca,'fontsize',12)
    text(.03,.9 ,[char(64+1) ')'],'units','normalized')

    subplot(3, 1, 2)
    plot(t, rel_phase_wrapped, 'k-', 'LineWidth', 2, 'DisplayName', 'Leader - Follower')
    set(gca,'YTick',[-pi/4 0 pi/4],'YtickLabel',[-pi/4 0 pi/4]./2/pi*360)
    ylim([-pi/4 pi/4])
    legend('hide')
    xlabel('Time, s')
    ylabel('\phi [degrees]')
    set(gca,'fontsize',12)
    text(.03,.9 ,[char(64+2) ')'],'units','normalized')

    subplot(3, 1, 3)
    plot(t, omega_m./2/pi, 'k-', 'LineWidth', 2, 'DisplayName', 'Leader')
    hold on
    plot(t, omega_s./2/pi, 'm--', 'LineWidth', 2, 'DisplayName', 'Follower')
    hold off
    xlabel('Time [s]')
    ylabel('\omega [rev/s]')
    legend('show','location','east')
    set(gca,'fontsize',12)
    ylim([min(omega_s./2/pi)*.99 max(omega_s./2/pi)*1.01])
    text(.03,.9 ,[char(64+3) ')'],'units','normalized')

    if save_fig == 1
        set(gcf,'color','w')
        set(gcf, 'PaperPosition', [0 0 5 6])
        set(gcf, 'InvertHardcopy', 'off')
        f = fullfile(pwd,['theta_phi_omega_dynamics_' char(datetime('now','TimeZone','local','Format','y-MM-d-HHmmss')) '.jpeg']);
        print('-djpeg','-r600',f)
    end
end


%% Relative phase
if plotting_flag == 1
    figure(24527)
    subplot(2, 2, 1)
    plot(t, rel_phase, 'b-', 'LineWidth', 2)
    xlabel('Time [s]')
    ylabel('\theta_m(t) - \theta_s(t)')
    grid on;

    subplot(2, 2, 2)
    plot(t, theta_m + tau_s(2), 'k-', 'LineWidth', 2, ...
        'DisplayName', 'Leader \theta_m(t+\tau)')
    hold on
    plot(t, theta_s, 'm--', 'LineWidth', 2, 'DisplayName', 'Follower \theta_s(t)')
    hold off
    xlabel('Time [s]')
    ylabel('Verification: \theta_s(t) vs \theta_m(t+\tau)')
    legend('show')
    grid on

    subplot(2, 2, 3)
    plot(t, coupling_term, 'b-', 'LineWidth', 2)
    xlabel('Time [s]')
    ylabel('Coupling Term')
    grid on

    subplot(2, 2, 4)
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 0.5); hold on
    l(1) = plot(cos(theta_m(1)), sin(theta_m(1)), 'k^', 'MarkerSize', 8, 'LineWidth', 2);
    l(2) = plot(cos(theta_s(1)), sin(theta_s(1)), 'm^', 'MarkerSize', 8, 'LineWidth', 2);
    l(3) = plot(cos(theta_m(end)), sin(theta_m(end)), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
    l(4) = plot(cos(theta_s(end)), sin(theta_s(end)), 'ms', 'MarkerSize', 8, 'LineWidth', 2);
    hold off
    axis equal; grid on
    xlabel('cos(\phi)'); ylabel('sin(\phi)')
    legend(l,'Leader Start', 'Follower Start', 'Leader End', 'Follower End', 'location', 'southeast')

    set(gcf,'color','w')
    if save_fig == 1
        f = fullfile(pwd,['parameter_drift_dynamics_' char(datetime('now','TimeZone','local','Format','y-MM-d-HHmmss')) '.jpeg']);
        print('-djpeg','-r300',f)
    end
end

end

%% Run a single simulation
%  This returns PLV, mean relative phase, and (optionally) full traces
function [t_out, TH_L, TH_F, OMF, coupling_term] = ...
    run_sim(Omega_F, Theta_m0, Theta_s0, k1, k2, eps1, eps2, ...
    mu, alpha, tau_L, tau_F1, tau_F2, sigma, ...
    dt, nSteps)

% Delay samples
tauL_s  = max(1, round(tau_L  / dt));
tauF1_s = max(1, round(tau_F1 / dt));
tauF2_s = max(1, round(tau_F2 / dt));
buf_len = max([tauL_s, tauF1_s, tauF2_s]) + 2;

% Circular buffers for theta_L and theta_F histories
buf_L   = zeros(1, buf_len);
buf_F   = zeros(1, buf_len);
buf_ptr = 1;

% Initial state: [theta_L, theta_F, omega_F]
state      = [Theta_m0; Theta_s0; Omega_F];
buf_L(:)   = state(1);
buf_F(:)   = state(2);

% Outputs
counter  = 0;
t_out = (0:nSteps-1)*dt;
noise = randn(nSteps+3, 2)*pi*sigma;
TH_L  = zeros(nSteps, 1);
TH_F  = zeros(nSteps, 1);
OMF   = zeros(nSteps, 1);
coupling_term = zeros(nSteps, 1);

for k = 1:nSteps
    counter = counter + 1;

    TH_L(k) = state(1);
    TH_F(k) = state(2);
    OMF(k)  = state(3);
    coupling_term(k) = nan;

    % Read delayed values from buffers
    ptr_L   = mod(buf_ptr - tauL_s  - 1, buf_len) + 1;
    ptr_F1  = mod(buf_ptr - tauF1_s - 1, buf_len) + 1;
    ptr_F2  = mod(buf_ptr - tauF2_s - 1, buf_len) + 1;

    thL_del  = buf_L(ptr_L);
    thF_del1 = buf_F(ptr_F1);
    thF_del2 = buf_F(ptr_F2);

    % Write current state into buffers
    buf_L(buf_ptr) = state(1);
    buf_F(buf_ptr) = state(2);
    buf_ptr = mod(buf_ptr, buf_len) + 1;

    % RK4. Frozen delay history across stages
    k1s = dt * ode_rhs(state, thL_del, thF_del1, thF_del2, ...
        Omega_F, k1, k2, eps1, eps2, mu, alpha, noise(counter,:));
    k2s = dt * ode_rhs(state + k1s/2, thL_del, thF_del1, thF_del2, ...
        Omega_F, k1, k2, eps1, eps2, mu, alpha, noise(counter+1,:));
    k3s = dt * ode_rhs(state + k2s/2, thL_del, thF_del1, thF_del2, ...
        Omega_F, k1, k2, eps1, eps2, mu, alpha, noise(counter+2,:));
    k4s = dt * ode_rhs(state + k3s,   thL_del, thF_del1, thF_del2, ...
        Omega_F, k1, k2, eps1, eps2, mu, alpha, noise(counter+3,:));

    state = state + (k1s + 2*k2s + 2*k3s + k4s) / 6;
end
end


%% ODE right hand side
function ds = ode_rhs(state, thL_del, thF_del1, thF_del2, ...
    Omega_F, k1, k2, eps1, eps2, mu, alpha, noise)

th_L   = state(1);  %#ok<NASGU>
th_F   = state(2);  %#ok<NASGU>
omg_F  = state(3);

% Phase error terms (using delayed values)
phi1 = thL_del - thF_del1;
phi2 = thL_del - thF_del2;

% Leader: phase advances at current omega_F
dth_L  = omg_F + 0*noise(1,1);

% Follower: phase + coupling
coupling = k1 * sin(phi1) + k2 * sin(phi2);
dth_F  = omg_F + coupling + noise(1,2);

% Adaptive frequency (mu * domega_F/dt = ...)
domg_F = (1/mu) * ( alpha * (Omega_F - omg_F) ...
    + eps1 * sin(phi1) ...
    + eps2 * sin(phi2) );

ds = [dth_L; dth_F; domg_F];

end


%% Utility
function a = wrap_angle(a)
a = mod(a + pi, 2*pi) - pi;
end