function results = microEDM_nanoparticle_model(I_curr, V, freq, duty, varargin)
%% =========================================================================
%  NUMERICAL MODEL: NANOPARTICLE GENERATION BY MICRO-ELECTRO DISCHARGE
%  MACHINING (micro-EDM)
%
%  Implements the framework of:
%  Mitra, Muralidhara, Vasa & Singaperumal, "Investigation on particle
%  generation by micro-electro discharge machining", Proc. SPIE 7590,
%  75900J (2010).
%
%  Model structure (mirrors the paper's Sections 2.1-2.3 and Fig.6):
%    Thermal model of micro-EDM              (Sec. 2.1, eq 1-5)
%    Nucleation-coagulation + cooling-rate model, solved together in one
%    explicit time march                     (Sec. 2.2-2.3)
%
%  -------------------------------------------------------------------------
%  USAGE
%  -------------------------------------------------------------------------
%    results = microEDM_nanoparticle_model(I_curr, V, freq, duty)
%    results = microEDM_nanoparticle_model(I_curr, V, freq, duty, 'Name', Value, ...)
%
%  REQUIRED INPUTS
%    I_curr : average discharge current            [A]
%    V      : discharge voltage                    [V]
%    freq   : discharge frequency                  [Hz]
%    duty   : duty cycle, as a FRACTION (0 to 1), i.e. 0.3 = 30%   [-]
%             Pulse-on duration is derived internally as Ton = duty/freq.
%
%  OPTIONAL NAME-VALUE INPUTS (all have defaults matching the validated
%  steel-in-water case; override only what you need to change)
%    Material / workpiece properties:
%      'rho_p'      density                         [kg/m^3]   default 8000
%      'K_th'       thermal conductivity             [W/m-K]    default 12.5
%      'alpha_th'   thermal diffusivity              [m^2/s]    default 3.5e-6
%      'Tb'         boiling point                    [K]        default 3250
%      'Hlv'        latent heat of vaporization      [J/kg]     default 5.5e6
%      'Mmolar'     molar mass                       [kg/mol]   default 0.06617
%      'sigma'      surface tension (liquid metal)   [N/m]      default 1.8
%      'Cp_p'       specific heat capacity           [J/kg-K]   default 550
%      'Cfrac'      energy-transfer fraction         [-]        default 0.39
%    Dielectric fluid properties:
%      'Tf'         bulk fluid temperature           [K]        default 300
%      'mu_f'       dynamic viscosity                [Pa.s]     default 8.9e-4
%      'kf'         thermal conductivity             [W/m-K]    default 0.61
%      'Cp_f'       specific heat capacity           [J/kg-K]   default 4180
%      'rho_f'      density                          [kg/m^3]   default 997
%      'lambda_mfp' fluid mean free path              [m]        default 4e-10
%    Nodal (sectional) grid:
%      'd_min_nm'   smallest node diameter           [nm]       default 1
%      'd_max_nm'   largest node diameter            [nm]       default 200
%      'q_ratio'    node-to-node volume ratio        [-]        default 2
%    Simulation control:
%      'total_sim_time'  total simulated time [s]. Default = 1/freq (i.e.
%                         the inter-pulse period -- the model does not
%                         track interaction between successive pulses, so
%                         simulating beyond the next pulse's arrival has
%                         no physical basis; see the accompanying report
%                         notes on choosing this value).
%      'n_coarse_steps'  steps after the time-step has ramped up, default 3000
%      'n_ramp_steps'    steps used to ramp dt from fine to coarse, default 500
%      'Verbose'         print thermal-model & summary text, default true
%      'MakePlot'        create this run's own two figures, default false
%                         (leave false when batching multiple runs and
%                         overlaying them yourself -- see the companion
%                         example script compare_runs.m)
%
%  OUTPUT: a struct `results` with fields
%    .I_curr, .V, .freq, .duty, .Ton      -- the process parameters used
%    .R_um, .T_peak, .Ti                  -- thermal model results
%    .t          -- time vector [s]
%    .T          -- temperature history [K]
%    .d_avg      -- average particle diameter history [nm]
%    .d_star     -- critical diameter history [nm] (NaN once nucleation stops)
%    .label      -- a ready-made legend string for this run
%
%  -------------------------------------------------------------------------
%  MODELLING NOTES / ENGINEERING CORRECTIONS TO THE PUBLISHED EQUATIONS
%  -------------------------------------------------------------------------
%  The published paper is a scanned/OCR'd PDF and a few of its printed
%  equations are not directly usable as typeset. The following standard,
%  physically-consistent forms are used in their place so the code runs
%  and gives sensible numbers. This does not change the modelling
%  *framework*, only fixes typographical issues in specific formulas:
%
%   1) Eq.(13) Clausius-Clapeyron relation: the printed form
%      Ps = P0*exp[Hlv*(T-Tb)/(R*Tb)] is dimensionally inconsistent (the
%      exponent is not dimensionless) and overflows double precision for
%      the temperatures involved. The standard, dimensionally-consistent
%      Clausius-Clapeyron form is used instead:
%           Ps(T) = P0 * exp( Hlv*(T-Tb) / (Rgas*T*Tb) )
%
%   2) Eq.(1) heat-flux-radius correlation (Yeo, Kurnia & Tan, 2007):
%      current i in Amps and pulse-on time Ton in SECONDS are plugged
%      directly into R_um = 2040*i^0.43*Ton^0.44, and the numeric RESULT
%      is in MICROMETERS (the empirical constant absorbs the unit
%      conversion). This matches realistic micro-EDM crater sizes.
%
%   3) Eq.(14)-(15) nucleation rate: the paper's non-dimensional surface
%      tension term reuses the symbol used elsewhere for saturation
%      ratio, which is ambiguous/inconsistent. The standard classical
%      nucleation theory (CNT) form is used: phi = sigma*a1/(kb*T), where
%      a1 is the surface area of a single monomer (atom).
%
%   4) Eq.(7) mean free path of the dielectric fluid is not derived
%      anywhere in the paper; it is exposed here as an editable input
%      (lambda_mfp) with a reasonable default for water.
%
%   5) Reynolds number (used with eq.6, Nu ~ Re^2*Pr^2) is computed with
%      the standard definition Re = rho_f*v*d_p/mu_f, which the paper
%      does not spell out explicitly.
%
%   6) dt_min (the finest time step used at the start of the simulation)
%      is fixed in absolute terms (1 ps) rather than scaled from
%      total_sim_time. Scaling it caused numerical overshoot -- the
%      temperature would jump straight to the fluid temperature in a
%      single step -- whenever total_sim_time was set much larger than a
%      few microseconds.
%  =========================================================================

%% ---- 0. Parse optional Name-Value overrides ----------------------------
defaults = struct( ...
    'Cfrac',          0.39, ...
    'rho_p',          8000, ...
    'K_th',           12.5, ...
    'alpha_th',       3.5e-6, ...
    'Tb',             3250, ...
    'Hlv',            5.5e6, ...
    'Mmolar',         0.06617, ...
    'sigma',          1.8, ...
    'Cp_p',           550, ...
    'Tf',             300, ...
    'mu_f',           8.9e-4, ...
    'kf',             0.61, ...
    'Cp_f',           4180, ...
    'rho_f',          997, ...
    'lambda_mfp',     4e-10, ...
    'd_min_nm',       1, ...
    'd_max_nm',       200, ...
    'q_ratio',        2, ...
    'total_sim_time', [], ...      % empty => auto = 1/freq
    'n_coarse_steps', 3000, ...
    'n_ramp_steps',   500, ...
    'Verbose',        true, ...
    'MakePlot',       false);

opt = local_parse_namevalue(defaults, varargin);

Cfrac = opt.Cfrac; rho_p = opt.rho_p; K_th = opt.K_th; alpha_th = opt.alpha_th;
Tb = opt.Tb; Hlv = opt.Hlv; Mmolar = opt.Mmolar; sigma = opt.sigma; Cp_p = opt.Cp_p;
Tf = opt.Tf; mu_f = opt.mu_f; kf = opt.kf; Cp_f = opt.Cp_f; rho_f = opt.rho_f;
lambda_mfp = opt.lambda_mfp;
d_min_nm = opt.d_min_nm; d_max_nm = opt.d_max_nm; q_ratio = opt.q_ratio;
n_coarse_steps = opt.n_coarse_steps; n_ramp_steps = opt.n_ramp_steps;

% ---- Physical constants (fixed) -----------------------------------------
kb        = 1.380649e-23;  % Boltzmann constant           [J/K]
Na        = 6.02214076e23; % Avogadro's number             [1/mol]
P0        = 1.01325e5;     % Ambient / reference pressure  [Pa]
Rgas_univ = 8.314;         % Universal gas constant        [J/mol-K]

% ---- Process parameters --------------------------------------------------
Ton = duty / freq;                      % pulse-on duration [s], from duty cycle & frequency

if isempty(opt.total_sim_time)
    total_sim_time = 1/freq;            % default: one full inter-pulse period
else
    total_sim_time = opt.total_sim_time;
end

%% ---- SECTION 2: Derived parameters & nodal grid setup -------------------
Rgas_p = Rgas_univ / Mmolar;
m1     = Mmolar / Na;
v1     = m1 / rho_p;
a1     = (36*pi)^(1/3) * v1^(2/3);

v_min = (pi/6)*(d_min_nm*1e-9)^3;
v_max = (pi/6)*(d_max_nm*1e-9)^3;
n_nodes = ceil( log(v_max/v_min) / log(q_ratio) ) + 1;
v_node  = v_min * q_ratio.^(0:n_nodes-1)';
d_node  = (6*v_node/pi).^(1/3); %#ok<NASGU>

[idxK, idxKp1, wK, wKp1] = local_build_coagulation_operator(v_node);
G = local_build_collision_geometry(v_node);
const_beta = (3/(4*pi))^(1/6);

dt_min = 1e-12;                              % fixed fine time step [s]
dt_max = total_sim_time / n_coarse_steps;    % coarse time step, scales with sim length [s]
if dt_max < dt_min, dt_max = dt_min; end
growth = (dt_max/dt_min)^(1/n_ramp_steps);

%% ---- SECTION 3: Thermal model of micro-EDM ------------------------------
R_um = 2040 * I_curr^0.43 * Ton^0.44;        % [micrometers]
R    = R_um * 1e-6;                          % heat-flux radius [m]

q_heat = Cfrac*V*I_curr / (pi*R^2);          % heat flux [W/m^2], eq.(2)

tau = alpha_th*Ton / R^2;                    % non-dimensional time, eq.(5)

theta_peak = 2*sqrt(tau/pi)*(1 - exp(-1/(4*tau))) + erfc(1/(2*sqrt(tau))); % eq.(3)

T_peak = theta_peak * (Cfrac*V*I_curr) / (K_th*pi*R);  % eq.(4) rearranged

Ti = 0.5*(T_peak + Tb);

if opt.Verbose
    fprintf('--- Thermal model results (I=%gA, V=%gV, f=%gHz, duty=%g) ---\n', I_curr, V, freq, duty);
    fprintf('Pulse-on duration Ton   = %.3f us\n', Ton*1e6);
    fprintf('Heat flux radius R      = %.3f um\n', R_um);
    fprintf('Heat flux q             = %.3e W/m^2\n', q_heat);
    fprintf('Peak temperature T_peak = %.1f K\n', T_peak);
    fprintf('Initial vapor temp Ti   = %.1f K\n\n', Ti);
end

%% ---- SECTIONS 4-5: Nucleation-coagulation + cooling-rate model ----------
T  = Ti;
Ps_i = local_clausius_clapeyron(Ti, P0, Hlv, Rgas_p, Tb);
N1 = Ps_i / (kb*Ti);
Nk = zeros(n_nodes,1);

t  = 0;
dt = dt_min;

maxsteps   = n_ramp_steps + n_coarse_steps + 10;
t_hist     = zeros(maxsteps,1);
Tt_hist    = zeros(maxsteps,1);
davg_hist  = zeros(maxsteps,1);
dstar_hist = NaN(maxsteps,1);
step = 0;

while t < total_sim_time && step < maxsteps
    step = step + 1;

    % ---- Saturation state & critical diameter (eqs. 12-13) ----
    Ps = local_clausius_clapeyron(T, P0, Hlv, Rgas_p, Tb);
    Ns = Ps/(kb*T);
    S  = N1/Ns;

    if S > 1
        d_star = 4*sigma*v1 / (kb*T*log(S));
        v_star = (pi/6)*d_star^3;

        phi = sigma*a1/(kb*T);
        Jk  = Ns*N1*v1*sqrt(2*sigma/(pi*m1)) * exp( phi - 4*phi^3/(27*log(S)^2) );

        knuc = find(v_node >= v_star, 1, 'first');
        if isempty(knuc), knuc = n_nodes; end
        eps_k = min(v_star / v_node(knuc), 1);
    else
        d_star = NaN; v_star = 0; Jk = 0; knuc = []; eps_k = 0;
    end

    % ---- Coagulation (eqs. 16, 20-21) ----
    beta  = const_beta * sqrt(6*kb*T/rho_p) * G;
    Mmat  = beta .* (Nk*Nk');
    Birth = accumarray(idxK,   0.5*wK  .*Mmat(:), [n_nodes,1]) + ...
            accumarray(idxKp1, 0.5*wKp1.*Mmat(:), [n_nodes,1]);
    Death = Nk .* (beta*Nk);
    dNk_coag = Birth - Death;

    dNk_nuc = zeros(n_nodes,1);
    if ~isempty(knuc)
        dNk_nuc(knuc) = Jk*eps_k;
    end

    % ---- Update particle-size distribution ----
    Nk = max(Nk + dt*(dNk_nuc + dNk_coag), 0);
    N1 = max(N1 - dt*Jk*(v_star/v1), 0);

    % ---- Average particle diameter ----
    if sum(Nk) > 0
        v_avg = sum(Nk.*v_node)/sum(Nk);
        d_avg = (6*v_avg/pi)^(1/3);
    else
        d_avg = 0;
    end

    % ---- Cooling-rate model (eqs. 6-11) ----
    if d_avg > 0
        d_p = d_avg;
    elseif ~isnan(d_star)
        d_p = d_star;
    else
        d_p = (6*v1/pi)^(1/3);
    end
    Pr    = mu_f*Cp_f/kf;
    Do    = kb*T / (3*pi*mu_f*d_p);
    v_vel = Do/lambda_mfp;
    Re    = rho_f*v_vel*d_p/mu_f;
    Nu    = Re^2 * Pr^2;
    h     = Nu*kf/d_p;
    coolrate = 6*h/(rho_p*Cp_p*d_p) * (T - Tf);
    T = max(T - coolrate*dt, Tf);

    % ---- store history ----
    t_hist(step)     = t;
    Tt_hist(step)    = T;
    davg_hist(step)  = d_avg*1e9;
    dstar_hist(step) = d_star*1e9;

    t  = t + dt;
    dt = min(dt*growth, dt_max);
end

t_hist     = t_hist(1:step);
Tt_hist    = Tt_hist(1:step);
davg_hist  = davg_hist(1:step);
dstar_hist = dstar_hist(1:step);

if opt.Verbose
    fprintf('--- Simulation summary ---\n');
    fprintf('Final temperature reached : %.1f K\n', Tt_hist(end));
    fprintf('Final average diameter    : %.2f nm\n', davg_hist(end));
    fprintf('Minimum critical diameter : %.3f nm\n\n', min(dstar_hist));
end

%% ---- Package outputs -----------------------------------------------------
results.I_curr = I_curr;
results.V      = V;
results.freq   = freq;
results.duty   = duty;
results.Ton    = Ton;
results.R_um   = R_um;
results.T_peak = T_peak;
results.Ti     = Ti;
results.t      = t_hist;
results.T      = Tt_hist;
results.d_avg  = davg_hist;
results.d_star = dstar_hist;
results.label  = sprintf('I=%gA, V=%gV, f=%gHz, duty=%g', I_curr, V, freq, duty);

%% ---- Optional standalone plots for a single run --------------------------
if opt.MakePlot
    figure('Name','Average particle diameter vs time');
    plot(results.t*1e6, results.d_avg, 'b-', 'LineWidth', 1.8);
    xlabel('Time (\mus)'); ylabel('Average particle diameter (nm)');
    title(['Average nanoparticle diameter vs time  (' results.label ')']);
    grid on;

    figure('Name','Critical diameter vs time');
    plot(results.t*1e9, results.d_star, 'r-', 'LineWidth', 1.8);
    xlabel('Time (ns)'); ylabel('Critical diameter, d^* (nm)');
    title(['Critical diameter vs time  (' results.label ')']);
    grid on;
end

end % ======================= end of main function =========================


%% =========================================================================
%% LOCAL (SUB) FUNCTIONS
%% =========================================================================

function opt = local_parse_namevalue(defaults, args)
% Overrides fields of `defaults` with Name-Value pairs given in `args`,
% erroring if an unrecognized name is passed (helps catch typos).
    opt = defaults;
    if mod(numel(args),2) ~= 0
        error('microEDM_nanoparticle_model:badArgs', ...
            'Optional arguments must be given as Name-Value pairs.');
    end
    fn = fieldnames(defaults);
    for k = 1:2:numel(args)
        name = args{k};
        if ~ischar(name)
            error('microEDM_nanoparticle_model:badArgs', ...
                'Optional argument names must be strings.');
        end
        match = fn(strcmpi(fn, name));
        if isempty(match)
            error('microEDM_nanoparticle_model:unknownOption', ...
                'Unknown optional parameter "%s".', name);
        end
        opt.(match{1}) = args{k+1};
    end
end

function Ps = local_clausius_clapeyron(T, P0, Hlv, Rgas, Tb)
% Dimensionally-consistent Clausius-Clapeyron relation (corrected form of
% the paper's eq.13 -- see Modelling Note 1 at the top of this file).
    exponent = Hlv*(T - Tb) / (Rgas*T*Tb);
    exponent = min(exponent, 700);   % overflow guard (exp(700) ~ 1e304)
    Ps = P0*exp(exponent);
end

function [idxK, idxKp1, wK, wKp1] = local_build_coagulation_operator(v_node)
% Fixed-pivot (Hounslow) redistribution indices/weights used to split the
% volume of a coagulated pair (v_i+v_j) between the two bracketing nodes,
% per the paper's eq.(21).
    n = numel(v_node);
    idxK   = zeros(n*n,1);
    idxKp1 = zeros(n*n,1);
    wK     = zeros(n*n,1);
    wKp1   = zeros(n*n,1);
    c = 0;
    for i = 1:n
        for j = 1:n
            c = c + 1;
            vij = v_node(i) + v_node(j);
            if vij >= v_node(end)
                idxK(c) = n; idxKp1(c) = n; wK(c) = 1; wKp1(c) = 0;
            else
                k = find(v_node <= vij, 1, 'last');
                if k == n
                    idxK(c) = n; idxKp1(c) = n; wK(c) = 1; wKp1(c) = 0;
                else
                    kp1 = k + 1;
                    wk_ = (v_node(kp1) - vij) / (v_node(kp1) - v_node(k));
                    idxK(c)   = k;
                    idxKp1(c) = kp1;
                    wK(c)     = wk_;
                    wKp1(c)   = 1 - wk_;
                end
            end
        end
    end
end

function G = local_build_collision_geometry(v_node)
% Purely volume-dependent part of the Brownian collision frequency,
% eq.(16): beta(i,j) = const * sqrt(T) * G(i,j).
    n = numel(v_node);
    G = zeros(n,n);
    for i = 1:n
        for j = 1:n
            G(i,j) = sqrt(1/v_node(i) + 1/v_node(j)) * ...
                     (v_node(i)^(1/3) + v_node(j)^(1/3))^2;
        end
    end
end