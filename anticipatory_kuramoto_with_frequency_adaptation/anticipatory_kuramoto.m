function [mean_rel_phase_consistency,mean_rel_phase,rel_phase_error,tau,omega_delta,omega_s] = anticipatory_kuramoto(varargin)
%ANTICIPATORY_KURAMOTO  Anticipatory Kuramoto model with delayed feedback and frequency adaptation
%   [mean_rel_phase_consistency,mean_rel_phase,rel_phase_error,tau,omega_delta,omega_s] = anticipatory_kuramoto(varargin)
%   Arguments:
%      [K_1 and K_2]
%      [epsilon_1 epsilon_2]
%      [tau_1 tau_2]
%      tau_leader
%      relative_phase_initial
%      empty
%      show figs flag
%      print figs flag

%   Dobri Dotov, UNO, 2026


%% Parameters
dt = 0.001;                  % Time step
T = 6;                       % Total time
t = 0:dt:T;                  % Time vector
n = length(t);


% Natural frequencies
omega_m = zeros(n,1);      % Leader freq
omega_s = zeros(n,1);      % Follower freq
omega_m(1) = 2*2*pi;       % Leader freq 0
omega_s(1) = 2*2*pi;       % Follower freq 0
coupling_term = zeros(n,1);


%% Initialization
phi_m = zeros(1, n);    % Leader phase
phi_s = zeros(1, n);    % Follower phase


if numel(varargin) > 0
    K = varargin{1}; % Phase coupling strength
    epsilon = varargin{2}; % Frequency adaptation strength
    tau = varargin{3}; % Delay phase (anticipation horizon)
    tau_m = varargin{4};
    % Initial rel phase
    phi_m(1) = 0;
    phi_s(1) = -varargin{5};
    % Viz
    plotting_flag = varargin{7};
    save_fig = varargin{8};
else % Default parameters
    K = [1 1];
    epsilon = [1 1];
    tau = [.2*2*pi .2*2*pi];
    tau_m = 2*pi*.2;
    phi_m(1) = 0;
    phi_s(1) = 0;
    plotting_flag = 1;
    save_fig = 0;
end

% Delay parameters
delay_m = round(tau_m/omega_s(1)/dt);
delay_s(1) = round(tau(1)/omega_s(1)/dt);
delay_s(2) = round(tau(2)/omega_s(1)/dt);
mu_timescale = 10; % freq adaptation timescale
alpha = .0; % ignore for now
sigma = .0;


%% Euler loop for the Kuramoto coupling with delayed feedback and freq adaptation
for i = 2:n
    phi_m(i) = phi_m(i-1) + dt*omega_m(i-1); % Update master's simple rotation

    on = (i - max(delay_m,max(delay_s)) - 1)>0; % Wait for the delay period before enabling the coupling
    if on
        coupling_term(i) = K(1)*sin(phi_m(max(1,i - delay_m - 1)) - phi_s(i - delay_s(1) - 1)) + K(2)*sin(phi_m(max(1,i - delay_m - 1)) - phi_s(max(1,i - delay_s(2) - 1)));
    else
        coupling_term(i) = 0;
    end
    phi_s(i) = phi_s(i-1) + dt*(omega_s(i-1) + coupling_term(i));

    if on
        coupling_term(i) = epsilon(1)*sin(phi_m(max(1,i - delay_m - 1)) - phi_s(i - delay_s(1) - 1)) + epsilon(2)*sin(phi_m(max(1,i - delay_m - 1)) - phi_s(max(1,i - delay_s(2) - 1)));
    else
        coupling_term(i) = 0;
    end
    omega_s(i) = omega_s(i-1) + mu_timescale^-1*dt*(alpha*(omega_s(1) - omega_s(i-1)) + coupling_term(i) + sigma*randn*2*pi);
    omega_m(i) = omega_s(i);
end


%% Visualize
% Calculate synchronization measures
% Wrapped phase difference for plotting
% Phase difference
rel_phase_wrapped = mod(phi_m - phi_s + pi, 2*pi) - pi;
z = exp(1i * rel_phase_wrapped);
rel_phase = angle(z);
mean_rel_phase = angle(mean(z(round(numel(z)/2):end)));
mean_rel_phase_consistency = abs(mean(z(round(numel(z)/2):end)));
final_rel_phase = angle(mean(z(end-round(1/dt):end)));
omega_delta = (omega_s(end) - omega_s(1))/omega_s(1)*1e2;
rel_phase_error = mod(mean_rel_phase - tau_m + tau(1)*K(1)/sum(K) + tau(2)*K(2)/sum(K) + pi, 2*pi) - pi;


%% Display results
fprintf('External delay τ_m: %.2f rad\n', tau_m);
fprintf('Self-delay τ_1 and τ_2: %.2f and %.2f rad\n', tau);
fprintf('Phase coupling K_1 and K_2: %.2f and %.2f\n', K);
fprintf('Freq coupling ε_1 and ε_2: %.2f %.2f\n', epsilon);
fprintf('Final phase lag: %.6f rad\n', final_rel_phase);


%% Plotting
if plotting_flag == 1
    % Phase evolution
    figure(25421)
    subplot(3, 1, 1)
    plot(t, mod(phi_m,2*pi), 'k-', 'LineWidth', 2, 'DisplayName', 'Leader')
    hold on
    plot(t, mod(phi_s,2*pi), 'm-', 'LineWidth', 2, 'DisplayName', 'Follower')
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
    ylabel('\phi_m(t) - \phi_s(t)')
    grid on;

    subplot(2, 2, 2)
    plot(t, phi_m + tau(2), 'k-', 'LineWidth', 2, ...
        'DisplayName', 'Leader \phi_m(t+\tau)')
    hold on
    plot(t, phi_s, 'm--', 'LineWidth', 2, 'DisplayName', 'Follower \phi_s(t)')
    hold off
    xlabel('Time [s]')
    ylabel('Verification: \phi_s(t) vs \phi_m(t+\tau)')
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
    l(1) = plot(cos(phi_m(1)), sin(phi_m(1)), 'k^', 'MarkerSize', 8, 'LineWidth', 2);
    l(2) = plot(cos(phi_s(1)), sin(phi_s(1)), 'm^', 'MarkerSize', 8, 'LineWidth', 2);
    l(3) = plot(cos(phi_m(end)), sin(phi_m(end)), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
    l(4) = plot(cos(phi_s(end)), sin(phi_s(end)), 'ms', 'MarkerSize', 8, 'LineWidth', 2);
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
