%% bayes_gait_state_space_causal_policy_model.m
%
% A minimal predictive-processing / Bayesian causal-inference gait model.
%
% Purpose
% -------
% This script is a proof-of-principle simulation, not a fitted data model.
% It is designed to show how one relatively compact state-space model can
% generate two qualitative effects:
%
%   1. Amplified, non-delayed auditory feedback can reduce cadence
%      variability by increasing the precision of sensory evidence.
%
%   2. Delayed amplified auditory feedback can produce a nonlinear
%      cadence-delay relationship: small delays slow cadence, whereas larger
%      delays can speed cadence.
%
% Core theoretical idea
% ---------------------
% The observer/controller estimates the current step period T. On each step,
% the system combines:
%
%   A. an internal prediction/prior about the next step period,
%   B. undelayed natural multisensory evidence, and
%   C. experimentally amplified auditory evidence, which may be delayed.
%
% The key addition relative to simple cue integration is Bayesian causal
% gating. The delayed auditory signal is not always assumed to arise from the
% current step. Instead, the model computes a simple common-cause posterior:
%
%        L_common(d) = exp(-d^2 / (2*tau_c^2))
%
%                      pi_c * L_common(d)
%        p_common(d) = -------------------------------
%                      pi_c * L_common(d) + (1-pi_c)*eta
%
% Because this is a posterior probability rather than a raw temporal
% likelihood, p_common(0) is controlled by the common-cause prior and the
% separate-cause likelihood. With eta = 1, p_common(0) = pi_c.
%
% When delay is small, p_common is high, so the auditory signal is treated as
% self-generated sensory evidence and is integrated. When delay is large,
% p_common is low, so the auditory signal is downweighted as feedback from
% the current step. In that regime, the same sound can still act as an
% external cue / perturbation / attention-capturing event and trigger an
% "act-early" policy that shortens the next step period.

% Yi Gao, University of Nebraska Omaha, 2026

clear; close all
rng(1, 'twister');

%% Natural baseline auditory gain.
% This is intentionally > 0. The baseline condition is not a literal 
% no-sound condition; it represents ordinary walking where natural 
% footstep sounds are present but low-gain and in the background.
sim.ampNatural = 0.20;

%% Natural baseline condition
%
% This is the reference condition. It contains:
%   - internal prediction,
%   - undelayed natural multisensory evidence,
%   - weak natural auditory evidence with zero experimental delay.
%
% It is not labeled "no auditory" because the model assumes some baseline
% natural sound is always available.
theta = defaultTheta();
sim.nTrials = 80;      % repeated simulated trials per condition
sim.nSteps  = 1200;    % steps per simulated trial
sim.burnIn  = 200;     % discard early transient steps
baseline = runCondition(0, sim.ampNatural, theta, sim);
baseCad = baseline.cadMean;
baseCV  = baseline.cadCV;
baseIntervalCV = baseline.intervalCV;

%%
switch 1
    case 0
        selected_parameters_or_full_space_flag = 'selected'; % Quick demo.
    case 1
        selected_parameters_or_full_space_flag = 'full'; % The full space can take time to complete.
end

switch selected_parameters_or_full_space_flag
    case 'selected'
        % Representative experimental amplified-auditory conditions, not
        % the full space
        %
        % delayFracs are expressed as fractions of the current step period. A value
        % of 0.125 means the auditory feedback is delayed by 12.5% of the step
        % period. This matches the idea that auditory delay is relative to gait
        % timing, not an absolute perceptual event.
        delayFracs = [0, 0.125, 0.25];

        % amp = experimentally controlled auditory gain.
        % In this example, amp = 0.5 is moderate amplification and amp = 1.0 is high
        % amplification. Larger amp increases auditory precision.
        amps = [0.5, 1.0];

        meanCad = zeros(numel(amps), numel(delayFracs));
        meanCV  = zeros(numel(amps), numel(delayFracs));
        meanIntervalCV = zeros(numel(amps), numel(delayFracs));
        meanPCommon = zeros(numel(amps), numel(delayFracs));
        meanWAud = zeros(numel(amps), numel(delayFracs));
        meanPolicy = zeros(numel(amps), numel(delayFracs));
        for ia = 1:numel(amps)
            for id = 1:numel(delayFracs)
                out = runCondition(delayFracs(id), amps(ia), theta, sim);

                meanCad(ia,id) = out.cadMean;
                meanCV(ia,id)  = out.cadCV;
                meanIntervalCV(ia,id) = out.intervalCV;
                meanPCommon(ia,id) = out.pCommonMean;
                meanWAud(ia,id) = out.wAudMean;
                meanPolicy(ia,id) = out.policyMean;
            end
        end

        % Express outcomes relative to natural baseline.
        dCad = (meanCad - baseCad) ./ baseCad * 100;
        dCV  = (meanCV  - baseCV)  ./ baseCV  * 100;
        dIntervalCV = (meanIntervalCV - baseIntervalCV) ./ baseIntervalCV * 100;

        % Print compact numerical summary
        fprintf('\nNatural baseline: cadence = %.2f steps/min, CV = %.3f%%\n', ...
            baseCad, baseCV);
        fprintf('Natural baseline: interval CV = %.3f%%\n', baseIntervalCV);

        fprintf('\nCondition summaries relative to natural baseline:\n');
        fprintf('%12s','Amp','DelayFrac','dCad(%)','dCadCV(%)','dIntCV(%)','pCommon','wAud','policy(ms)');
        fprintf('\n')
        for ia = 1:numel(amps)
            for id = 1:numel(delayFracs)
                fprintf('%12.2f', ...
                    amps(ia), delayFracs(id), dCad(ia,id), dCV(ia,id), dIntervalCV(ia,id), ...
                    meanPCommon(ia,id), meanWAud(ia,id), meanPolicy(ia,id) * 1000);
                fprintf('\n')
            end
        end

        %% Plot 1: Cadence change
        figure('Color','w');
        plot(delayFracs, dCad(1,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
        plot(delayFracs, dCad(2,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        yline(0, '--', 'LineWidth', 1);
        xlabel('Auditory delay (fraction of step period)');
        ylabel('\Delta cadence (%) vs natural baseline');
        legend('Amp = 0.5', 'Amp = 1.0', 'Location', 'Best');
        title('Cadence: small-delay slowing and larger-delay speed-up');
        grid on;

        %% Plot 2: Cadence variability change
        figure('Color','w');
        plot(delayFracs, dCV(1,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
        plot(delayFracs, dCV(2,:), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        yline(0, '--', 'LineWidth', 1);
        xlabel('Auditory delay (fraction of step period)');
        ylabel('\Delta CV(cadence) (%) vs natural baseline');
        legend('Amp = 0.5', 'Amp = 1.0', 'Location', 'Best');
        title('Cadence variability: precision benefit strongest at zero delay');
        grid on;

        %% Plot 3: Latent model quantities
        % it shows why the model behaves as it does:
        %
        %   - p_common decreases with delay.
        %   - auditory feedback weight decreases as delay increases.
        %   - policy term becomes more negative at larger delays.
        %
        % A negative policy term shortens step period and therefore increases
        % cadence.
        figure('Color','w');
        tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

        nexttile;
        plot(delayFracs, meanPCommon(1,:), '-o', 'LineWidth', 1.8); hold on;
        plot(delayFracs, meanPCommon(2,:), '-o', 'LineWidth', 1.8);
        ylim([0 1.05]);
        xlabel('Delay fraction');
        ylabel('p(common cause)');
        title('Causal attribution');
        grid on;

        nexttile;
        plot(delayFracs, meanWAud(1,:), '-o', 'LineWidth', 1.8); hold on;
        plot(delayFracs, meanWAud(2,:), '-o', 'LineWidth', 1.8);
        ylim([0 1.05]);
        xlabel('Delay fraction');
        ylabel('Effective auditory weight');
        title('Auditory feedback weight');
        grid on;

        nexttile;
        plot(delayFracs, meanPolicy(1,:) * 1000, '-o', 'LineWidth', 1.8); hold on;
        plot(delayFracs, meanPolicy(2,:) * 1000, '-o', 'LineWidth', 1.8);
        yline(0, '--', 'LineWidth', 1);
        xlabel('Delay fraction');
        ylabel('Policy shift (ms)');
        title('Late-penalty policy');
        grid on;
        legend('Amp = 0.5', 'Amp = 1.0', 'Location', 'Best');
    case 'full'
        % Sample across a space of two parameter models
        % The grid scans auditory delay
        % and auditory amplification, then shows how the same fixed model
        % architecture maps those two experimental manipulations onto:
        %   - mean cadence, expressed as percent change from natural baseline
        %   - cadence variability, expressed as percent change in CV(cadence)
        sim.delayGrid = linspace(0, 1, 40);
        sim.ampGrid = linspace(0.20, 1.20, 40);
        spaceOut = runParameterSpace(sim.delayGrid, sim.ampGrid, theta, ...
            sim, baseCad, baseCV, baseIntervalCV);

        % creates the figure which is figure 7 in the manuscript.
        fig_param_space(sim, spaceOut)
end


%%
function theta = defaultTheta()
% defaultTheta
%
% Returns all model parameters. The values below are not fit to data. They
% are chosen to produce realistic-looking qualitative behavior while keeping
% the model transparent.

% Internal preferred step period.
% T0 = 0.5 seconds corresponds to 120 steps/min.
theta.T0 = 0.50;

% Initial uncertainty over the internal estimate of step period.
theta.P0 = (0.010)^2;

% Process noise: uncertainty added at each prediction step.
% This captures slow internal drift / unmodeled fluctuations in gait timing.
theta.Q = (0.002)^2;

% Undelayed natural multisensory evidence.
% This cue summarizes proprioception, haptic feedback, vision, vestibular
% input, and weak natural auditory input that are not experimentally delayed.
theta.sigmaNatural = 0.020;

% Auditory evidence noise.
% sigmaAud decreases as amplitude increases. Delay is not included directly
% in auditory variance here; instead, delay affects causal attribution
% through p_common and affects the auditory observation through audBias.
%
% Earlier versions used sigmaAud = sigmaAud0 / amp + delayNoise. That mapping
% can make the high-amplitude condition unrealistically precise and therefore
% overpredict the observed variability reduction. The saturating mapping
% below keeps amplitude beneficial, but bounded:
%
%   sigmaAud = floor + gain / (1 + ampPrecisionGain * amp)
%
% where "floor" is the irreducible auditory timing noise.
theta.sigmaAudFloor = 0.006;
theta.sigmaAudGain = 0.020;
theta.ampPrecisionGain = 1.5;

% Delay bias in auditory evidence.
% When the delayed sound is still treated as self-generated feedback, it pulls
% the estimated step period upward. A larger period means slower cadence.
theta.biasDelay = 0.25;

% Causal-attribution parameters.
%
% The common-cause likelihood is Gaussian-shaped in auditory delay:
%
%   L_common(d) = exp(-d^2 / (2*tau_c^2))
%
% This likelihood is converted into a normalized common-cause posterior:
%
%                      pi_c * L_common(d)
%   p_common(d) = -------------------------------
%                 pi_c * L_common(d) + (1-pi_c)*eta
%
% piCommonPrior is the prior probability that amplified footstep sound and
% the current step belong to the same self-feedback event before considering
% the measured delay. separateCauseLikelihood is the relative likelihood of
% the separate-cause model across the tested delay range. Setting eta = 1
% treats the separate-cause model as a broad, approximately flat temporal
% alternative.
% With eta = 1, p_common(0) equals piCommonPrior; increase piCommonPrior if
% zero-delay feedback should be treated as almost certainly self-generated.
%
% Smaller tauCommon means the model stops treating delayed sound as
% self-generated feedback at smaller delays. With the defaults below, the
% common-cause posterior is high at zero delay, remains moderate at 12.5% of
% the step period, and is strongly reduced at 25%.
theta.piCommonPrior = 0.90;
theta.tauCommon = 0.095;
theta.separateCauseLikelihood = 1.00;

% Late-action penalty strength.
% This is the main modeling hypothesis beyond standard cue integration. When
% delayed auditory feedback is unlikely to be feedback from the current step,
% the controller is assumed to treat being late as more costly than being
% early. This creates a small negative shift in chosen step period.
%
% lambdaLate = latePenaltyStrength * delayFrac * (1 - p_common)
%
% If lambdaLate = 0, early and late timing errors have equal cost and the
% optimal action is the posterior mean. If lambdaLate > 0, the optimal action
% is shifted earlier, with the magnitude scaled by posterior uncertainty.
theta.latePenaltyStrength = 8.0;

% Motor execution noise.
% Execution noise increases with posterior uncertainty.
theta.motor0 = 0.0015;
theta.motorScale = 0.45;

% Plausible bounds on step period.
theta.Tmin = 0.32;
theta.Tmax = 0.75;

% Model switches useful for ablation tests.
theta.useCausalGate = true;
theta.usePolicy = true;
end


function cond = runCondition(delayFrac, amp, theta, sim)
% Repeats simulateTrial several times and returns averaged summary measures.

cadMeanAll = zeros(sim.nTrials,1);
cadCVAll = zeros(sim.nTrials,1);
intervalCVAll = zeros(sim.nTrials,1);
pCommonAll = zeros(sim.nTrials,1);
wAudAll = zeros(sim.nTrials,1);
policyAll = zeros(sim.nTrials,1);

for t = 1:sim.nTrials
    trial = simulateTrial(delayFrac, amp, theta, sim);

    cadMeanAll(t) = trial.cadMean;
    cadCVAll(t) = trial.cadCV;
    intervalCVAll(t) = trial.intervalCV;
    pCommonAll(t) = trial.pCommonMean;
    wAudAll(t) = trial.wAudMean;
    policyAll(t) = trial.policyMean;
end

cond.cadMean = mean(cadMeanAll);
cond.cadCV = mean(cadCVAll);
cond.intervalCV = mean(intervalCVAll);
cond.pCommonMean = mean(pCommonAll);
cond.wAudMean = mean(wAudAll);
cond.policyMean = mean(policyAll);
end


function spaceOut = runParameterSpace(delayGrid, ampGrid, theta, sim, ...
    baseCad, baseCadCV, baseIntervalCV)
% Evaluates the fixed toy model over a two-dimensional delay-by-amplitude
% grid. This is a sensitivity visualization, not a fitted parameter search.

nAmp = numel(ampGrid);
nDelay = numel(delayGrid);

cadMean = zeros(nAmp, nDelay);
cadCV = zeros(nAmp, nDelay);
intervalCV = zeros(nAmp, nDelay);

totalIterations = nAmp * nDelay;
iteration = 0;


for ia = 1:nAmp
    for id = 1:nDelay
        out = runCondition(delayGrid(id), ampGrid(ia), theta, sim);
        cadMean(ia,id) = out.cadMean;
        cadCV(ia,id) = out.cadCV;
        intervalCV(ia,id) = out.intervalCV;

        iteration = iteration + 1;

        if mod(iteration,10) == 0
            fprintf('Completed %d / %d iterations (%.1f%%)\n', ...
                iteration, ...
                totalIterations, ...
                100 * iteration / totalIterations);
        end

    end
end

spaceOut.cadMean = cadMean;
spaceOut.cadCV = cadCV;
spaceOut.intervalCV = intervalCV;
spaceOut.dCad = (cadMean - baseCad) ./ baseCad * 100;
spaceOut.dCadCV = (cadCV - baseCadCV) ./ baseCadCV * 100;
spaceOut.dIntervalCV = (intervalCV - baseIntervalCV) ./ baseIntervalCV * 100;
end


function out = simulateTrial(delayFrac, amp, theta, sim)
% simulateTrial
%
% Simulates one walking trial using a scalar state-space model.
%
% State:
%   T_hat = internal estimate of step period.
%   P     = uncertainty/variance of that estimate.
%
% Observations:
%   z_nat = undelayed natural multisensory evidence.
%   z_aud = experimentally amplified auditory evidence.
%
% The auditory evidence is integrated only to the extent that it is reliable
% and causally attributed to the current step.

T_hat = theta.T0;
P = theta.P0;

T_exec = zeros(sim.nSteps,1);
pCommonTrace = zeros(sim.nSteps,1);
wAudTrace = zeros(sim.nSteps,1);
policyTrace = zeros(sim.nSteps,1);

for n = 1:sim.nSteps

    %% Step 1: Prediction
    %
    % Predict that the next step period will continue from the current
    % internal estimate. Add process noise to represent uncertainty that
    % accumulates before new sensory evidence arrives.

    T_prior = T_hat;
    P_prior = P + theta.Q;

    %% Step 2: Undelayed natural multisensory update

    % It should be interpreted broadly: proprioceptive, haptic, visual,
    % vestibular, and weak natural auditory information about the current
    % step period.

    z_nat = theta.T0 + theta.sigmaNatural * randn();
    R_nat = theta.sigmaNatural^2;

    [T_nat, P_nat] = gaussianUpdate(T_prior, P_prior, z_nat, R_nat);

    %% Step 3: Auditory evidence generation

    % Higher amplitude makes the auditory signal more precise. Larger delay
    % does not directly increase auditory variance in this simplified version;
    % it instead lowers causal attribution to the current step through the
    % p_common gate below.

    if amp > 0
        audBias = theta.biasDelay * delayFrac * theta.T0;
        sigmaAud = theta.sigmaAudFloor + ...
            theta.sigmaAudGain / (1 + theta.ampPrecisionGain * amp);
        R_aud = sigmaAud^2;

        z_aud = theta.T0 + audBias + sigmaAud * randn();

        %% Step 4: Bayesian causal gate
        %
        % p_common is the posterior probability that the auditory cue and the
        % current motor event share a common cause. The Gaussian-shaped delay
        % term is treated as the temporal likelihood of the common-cause
        % model, then normalized against a broad separate-cause alternative.
        % This is still a compact toy approximation, but it follows the same
        % model-comparison logic as Bayesian causal-inference accounts.

        if theta.useCausalGate
            p_common = commonCausePosterior(delayFrac, theta);
        else
            p_common = 1;
        end

        %% Step 5: Reliability-weighted auditory update
        %
        % First compute what the posterior would be if the auditory cue were
        % fully attributed to the current step.

        [T_fused, P_fused] = gaussianUpdate(T_nat, P_nat, z_aud, R_aud);

        % Then blend the fused estimate with the natural-evidence-only
        % estimate. This is a simple form of Bayesian causal model averaging:
        %
        %   common cause:    use fused natural + auditory estimate
        %   separate causes: use natural-evidence-only estimate

        T_post = p_common * T_fused + (1 - p_common) * T_nat;

        % Variance blending includes the between-model uncertainty term. This
        % is a standard mixture-variance calculation:
        %
        %   Var(X) = E[Var(X | model)] + Var(E[X | model])

        P_post = p_common * (P_fused + (T_fused - T_post)^2) + ...
            (1 - p_common) * (P_nat + (T_nat - T_post)^2);

        % Motor execution noise is scaled by within-model sensory
        % uncertainty, not by the full between-model causal ambiguity. The
        % causal ambiguity is still retained in P_post for the decision rule
        % below. This keeps high-delay conditions from becoming unrealistically
        % variable simply because the model entertains two causal explanations.
        P_exec = p_common * P_fused + (1 - p_common) * P_nat;

        % Diagnostic effective auditory weight.
        % This reflects the actual influence of the auditory observation on the
        % final model-averaged state estimate:
        K_aud = P_nat / (P_nat + R_aud);
        w_aud_effective = p_common * K_aud;
        %% Step 6: Asymmetric late-penalty action policy
        %
        % When the delayed auditory signal is unlikely to be feedback from the
        % current step, it can still affect action as an external cue or
        % perturbation. Instead of adding a hand-shaped delay^2 policy term,
        % this version uses a simple Bayesian decision rule.
        %
        % The posterior over step period is approximated as:
        %
        %   T ~ Normal(T_post, P_post)
        %
        % The chosen action a is the step period that minimizes expected
        % asymmetric quadratic loss:
        %
        %   L(a,T) = c_late  * (a - T)^2, if a > T  (late step)
        %          = c_early * (a - T)^2, if a <= T (early step)
        %
        % where c_late = 1 + lambdaLate and c_early = 1. The effective late
        % penalty is:
        %
        %   lambdaLate = latePenaltyStrength * delayFrac * (1 - p_common)
        %
        % Thus, no-delay and high-attribution conditions have no early-action
        % bias. High-delay and low-attribution conditions place more cost on
        % being late and shift the chosen action earlier. The returned
        % policyShift is negative when the optimal action is earlier than the
        % posterior mean.
        %
        % This is deliberately separated from Bayesian sensory fusion:
        %   - sensory fusion estimates what the current step period is;
        %   - policy changes how the next action is executed.
        %
        % That separation is important because high-delay speed-up is not a
        % straightforward prediction of reliability-weighted cue integration.

        if theta.usePolicy
            policyShift = asymmetricLatePenaltyShift(P_post, p_common, ...
                delayFrac, theta);
        else
            policyShift = 0;
        end
    else
        % If amp is zero, no experimentally amplified auditory cue is present.
        % The model still has natural multisensory evidence through z_nat.

        T_post = T_nat;
        P_post = P_nat;
        P_exec = P_nat;
        p_common = 0;
        w_aud_effective = 0;
        policyShift = 0;
    end

    %% Step 7: Motor execution
    %
    % The executed period is the posterior estimate plus the policy shift and
    % motor noise. Motor noise scales with P_exec, a within-model sensory
    % uncertainty term. The broader P_post, which also includes causal
    % ambiguity, is used for action selection above.

    sigmaExec = theta.motor0 + theta.motorScale * sqrt(P_exec);
    Tn = T_post + policyShift + sigmaExec * randn();

    % Keep periods in a plausible gait range.
    Tn = min(theta.Tmax, max(theta.Tmin, Tn));

    %% Step 8: Propagation
    %
    % The executed step becomes the starting point for the next internal
    % prediction. This makes the model a closed perception-action loop.

    T_exec(n) = Tn;
    pCommonTrace(n) = p_common;
    wAudTrace(n) = w_aud_effective;
    policyTrace(n) = policyShift;

    T_hat = Tn;
    P = P_post;
end

%% Summaries after burn-in

idx = (sim.burnIn + 1):sim.nSteps;
T_use = T_exec(idx);
cad = 60 ./ T_use;

out.cadMean = mean(cad);
out.cadCV = std(cad) / mean(cad) * 100;
out.intervalCV = std(T_use) / mean(T_use) * 100;
out.pCommonMean = mean(pCommonTrace(idx));
out.wAudMean = mean(wAudTrace(idx));
out.policyMean = mean(policyTrace(idx));
end

function policyShift = asymmetricLatePenaltyShift(P_post, p_common, delayFrac, theta)
% asymmetricLatePenaltyShift
%
% Computes the action shift implied by an asymmetric quadratic loss.
%
% Posterior belief:
%   T ~ Normal(T_post, P_post)
%
% Loss:
%   L(a,T) = c_late  * (a - T)^2, if a > T  (chosen period is too long)
%          = c_early * (a - T)^2, if a <= T (chosen period is too short)
%
% Because the loss is translation-invariant, the optimal action shift
% relative to T_post depends only on posterior standard deviation and on the
% cost ratio. The function returns:
%
%   policyShift = a_star - T_post
%
% A negative value means the chosen period is shorter than the posterior
% mean, corresponding to earlier stepping and increased cadence.
%
% The late penalty is intentionally tied to causal attribution:
%
%   lambdaLate = latePenaltyStrength * delayFrac * (1 - p_common)
%
% This means the policy does not activate when delay is zero or when auditory
% feedback is still strongly attributed to the current step.

sigmaPost = sqrt(max(P_post, 0));
lambdaLate = theta.latePenaltyStrength * delayFrac * (1 - p_common);

if sigmaPost == 0 || lambdaLate <= 0
    policyShift = 0;
    return;
end

cEarly = 1;
cLate = 1 + lambdaLate;

% Solve for y = (a_star - T_post) / sigmaPost. The first-order condition for
% the asymmetric quadratic loss is:
%
%   y * [cLate*Phi(y) + cEarly*(1-Phi(y))]
%       + phi(y) * (cLate - cEarly) = 0
%
% where Phi and phi are the standard normal CDF and PDF. We use bisection
% rather than fzero so that the script does not require toolboxes.
lo = -6;
hi = 6;
for iter = 1:60
    mid = (lo + hi) / 2;
    if latePenaltyDerivative(mid, cEarly, cLate) > 0
        hi = mid;
    else
        lo = mid;
    end
end

yStar = (lo + hi) / 2;
policyShift = sigmaPost * yStar;
end

function p_common = commonCausePosterior(delayFrac, theta)
% commonCausePosterior
%
% Computes a normalized common-cause posterior for the auditory feedback
% stream. This replaces the simpler unnormalized gate
% exp(-0.5*(delay/tau)^2).
%
% The common-cause model says that small auditory delays are more likely if
% the sound came from the current step:
%
%   L_common(d) = exp(-d^2 / (2*tau_c^2))
%
% The separate-cause model is represented as a broad temporal alternative
% with approximately constant relative likelihood eta:
%
%   L_separate(d) = eta
%
% Bayes rule gives:
%
%                      pi_c * L_common(d)
%   p_common(d) = -------------------------------
%                 pi_c * L_common(d) + (1-pi_c)*eta
%
% This is not meant to be a full fitted temporal-binding model. It is a
% transparent approximation that makes the causal-attribution assumption
% explicit and easy to vary in ablation/sensitivity checks.

tau = max(theta.tauCommon, eps);
pi_c = min(max(theta.piCommonPrior, eps), 1 - eps);
eta = max(theta.separateCauseLikelihood, eps);

commonLikelihood = exp(-0.5 * (delayFrac / tau)^2);
weightedCommon = pi_c * commonLikelihood;
weightedSeparate = (1 - pi_c) * eta;

p_common = weightedCommon / (weightedCommon + weightedSeparate);
end

function g = latePenaltyDerivative(y, cEarly, cLate)
% latePenaltyDerivative
%
% Standardized derivative of expected asymmetric quadratic loss.

Phi = 0.5 * (1 + erf(y / sqrt(2)));
phi = exp(-0.5 * y^2) / sqrt(2 * pi);
g = y * (cLate * Phi + cEarly * (1 - Phi)) + ...
    phi * (cLate - cEarly);
end

function [m_post, P_post] = gaussianUpdate(m_prior, P_prior, z, R)
% gaussianUpdate
%
% Scalar Kalman/Gaussian update.
%
% m_prior, P_prior:
%   prior mean and variance of the current step-period estimate.
%
% z, R:
%   sensory observation and observation variance.
%
% K:
%   Kalman gain. The gain is large when sensory evidence is more reliable
%   than the prior and small when sensory evidence is noisy.

K = P_prior / (P_prior + R);
m_post = m_prior + K * (z - m_prior);
P_post = (1 - K) * P_prior;
end


function fig_param_space(space,spaceOut)

print_fig = 0;

%% delta cadence
[X,Y] = meshgrid(space.delayGrid, space.ampGrid-space.ampNatural);
Z = spaceOut.dCad; % delta cadence

figure('Color','w','Position',[100 100 900 650]);

hold on

surf(X,Y,Z,...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','none');

% Contour lines on the surface
contour3(X,Y,Z,...
    12,...
    'k',...
    'LineWidth',0.8);

% Zero reference plane
mesh(X,Y,...
    zeros(size(Z)),...
    'EdgeColor',[0.35 0.35 0.35],...
    'FaceColor','none',...
    'LineWidth',0.5);


% Zero contour projected onto the reference plane
contour3(X,Y,Z,...
    [0 0],...
    'c',...
    'LineWidth',3);

colormap(cool(256));

% Symmetric color limits around zero
cadLim = max(abs(Z(:)));
clim([-cadLim cadLim]);

% Labels
xlabel('d','FontSize',16)
ylabel('A','FontSize',16)
zlabel('\Delta Cadence (%)','FontSize',16)

% Axis appearance
grid on
box on

set(gca,...
    'FontSize',13,...
    'LineWidth',1.2,...
    'Projection','perspective')

axis tight

hold off

view([324.02 26.38])

if print_fig == 1
    f = fullfile(pwd,['tempo_3d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])
end


%% delta cv
[X,Y] = meshgrid(space.delayGrid, space.ampGrid-space.ampNatural);
Z = spaceOut.dCadCV;

figure('Color','w','Position',[100 100 900 650]);

hold on

surf(X,Y,Z,...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','none');

% Contour lines on the surface
contour3(X,Y,Z,...
    12,...
    'k',...
    'LineWidth',0.8);

% Zero reference plane
mesh(X,Y,...
    zeros(size(Z)),...
    'EdgeColor',[0.35 0.35 0.35],...
    'FaceColor','none',...
    'LineWidth',0.5);

% Zero contour projected onto the reference plane
contour3(X,Y,Z,...
    [0 0],...
    'c',...
    'LineWidth',2);

% Custom cyan -> purple colormap
cmap = [
    0.10 0.85 0.90
    0.25 0.75 0.95
    0.45 0.60 0.95
    0.60 0.45 0.95
    0.75 0.25 0.95
    0.95 0.05 0.90];

colormap(interp1(...
    linspace(0,1,size(cmap,1)),...
    cmap,...
    linspace(0,1,256)));

% Symmetric color limits around zero
cadLim = max(abs(Z(:)));
clim([-cadLim cadLim]);

% Labels
xlabel('d','FontSize',16)
ylabel('A','FontSize',16)
zlabel('\Delta CV (%)','FontSize',16)

% Lighting
camlight headlight
lighting gouraud
material([0.8 0.2 0.1 10])

% Axis appearance
grid on
box on

set(gca,...
    'FontSize',13,...
    'LineWidth',1.2,...
    'Projection','perspective')

axis tight

hold off

view([302.3 41.4])

if print_fig == 1
    f = fullfile(pwd,['cv_3d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])
end

end