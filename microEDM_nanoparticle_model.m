%% =======================================================================
%  NUMERICAL MODEL: NANOPARTICLE GENERATION BY MICRO-ELECTRO DISCHARGE
%  MACHINING (micro-EDM)
%
%  Implements the framework of:
%  Mitra, Muralidhara, Vasa & Singaperumal, "Investigation on particle
%  generation by micro-electro discharge machining", Proc. SPIE 7590,
%  75900J (2010).
%
%  Model structure (mirrors the paper's Sections 2.1-2.3 and Fig.6):
%    Section 1 : INPUT PARAMETERS            -- edit this section only
%    Section 2 : Derived parameters & nodal (sectional) grid setup
%    Section 3 : Thermal model of micro-EDM              (Sec. 2.1, eq 1-5)
%    Section 4/5: Nucleation-coagulation + cooling-rate model, solved
%                 together in one explicit time march     (Sec. 2.2-2.3)
%    Section 6 : Results & plots
%
%  OUTPUTS:
%    Fig 1: Average nanoparticle diameter vs time (microseconds)
%    Fig 2: Critical diameter (d*) vs time (nanoseconds)
%
%  -----------------------------------------------------------------------
%  MODELLING NOTES / ENGINEERING CORRECTIONS TO THE PUBLISHED EQUATIONS
%  -----------------------------------------------------------------------
%  The published paper is a scanned/OCR'd PDF and a few of its printed
%  equations are not directly usable as typeset. The following standard,
%  physically-consistent forms were used in their place so the code runs
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
%  =========================================================================

clear; clc; close all;

%% =========================================================================
%% SECTION 1: INPUT PARAMETERS  --  EDIT THIS SECTION ONLY
%% =========================================================================

% ---- 1.1 Micro-EDM process parameters ----------------------------------
V        = 40;        % Discharge voltage                        [V]
I_curr   = 3;           % Average discharge current                 [A]
freq     = 6000;        % Frequency                 [Hz]
duty     = 0.3;         % Duty Cycle                 [%]
Ton      = duty/freq;       % Pulse-on duration                         [s]
Cfrac    = 0.39;        % Fraction of discharge energy transferred to workpiece [-]

% ---- 1.2 Workpiece / particle material properties (default: Aluminium) --
rho_p    = 8000;        % Density                                   [kg/m^3]
K_th     = 12.5;         % Thermal conductivity                      [W/m-K]
alpha_th = 3.5e-6;      % Thermal diffusivity                       [m^2/s]
Tb       = 3250;        % Boiling point                             [K]
Hlv      = 5.5e6;      % Latent heat of vaporization                [J/kg]
Mmolar   = 0.06617;    % Molar mass                                [kg/mol]
sigma    = 1.8;        % Surface tension of the liquid metal        [N/m]
Cp_p     = 550;         % Specific heat capacity                     [J/kg-K]

% ---- 1.3 Dielectric fluid properties (default: deionized water) --------
Tf         = 300;       % Bulk dielectric fluid temperature          [K]
mu_f       = 8.9e-4;    % Dynamic viscosity                          [Pa.s]
kf         = 0.61;       % Thermal conductivity                       [W/m-K]
Cp_f       = 4180;      % Specific heat capacity                     [J/kg-K]
rho_f      = 997;       % Density                                    [kg/m^3]
lambda_mfp = 4e-10;     % Mean free path of fluid molecules (order-of-magnitude estimate) [m]

% ---- 1.4 Nodal (sectional) method discretization ------------------------
d_min_nm = 1;           % smallest node particle diameter            [nm]
d_max_nm = 200;        % largest node particle diameter (1 micron)  [nm]
q_ratio  = 2;           % volume ratio between consecutive nodes (Hounslow grid)

% ---- 1.5 Simulation time control ----------------------------------------
total_sim_time = 1/freq;  % total physical time simulated              [s]
n_coarse_steps = 3000000;  % approx. number of steps after dt has ramped up
n_ramp_steps   = 500;   % number of steps used to ramp dt from fine to coarse

% ---- 1.6 Physical constants (do not normally need editing) -------------
kb        = 1.380649e-23;  % Boltzmann constant           [J/K]
Na        = 6.02214076e23; % Avogadro's number             [1/mol]
P0        = 1.01325e5;     % Ambient / reference pressure  [Pa]
Rgas_univ = 8.314;         % Universal gas constant        [J/mol-K]

%% =========================================================================
%% SECTION 2: DERIVED PARAMETERS & NODAL GRID SETUP  (automatic)
%% =========================================================================
Rgas_p = Rgas_univ / Mmolar;         % specific gas constant of the vapor [J/kg-K]
m1     = Mmolar / Na;                % mass of a single atom (monomer)    [kg]
v1     = m1 / rho_p;                 % monomer volume                      [m^3]
a1     = (36*pi)^(1/3) * v1^(2/3);   % monomer surface area (sphere)       [m^2]

v_min = (pi/6)*(d_min_nm*1e-9)^3;
v_max = (pi/6)*(d_max_nm*1e-9)^3;
n_nodes = ceil( log(v_max/v_min) / log(q_ratio) ) + 1;
v_node  = v_min * q_ratio.^(0:n_nodes-1)';      % [n_nodes x 1] node volumes   [m^3]
d_node  = (6*v_node/pi).^(1/3);                  % [n_nodes x 1] node diameters [m]

% ---- Precompute coagulation redistribution weights (fixed-pivot / Hounslow method, eq.21) ----
[idxK, idxKp1, wK, wKp1] = build_coagulation_operator(v_node);

% ---- Precompute purely geometric part of the collision-frequency matrix, eq.(16) ----
G = build_collision_geometry(v_node);
const_beta = (3/(4*pi))^(1/6);

% ---- Time-step ramp schedule (fine steps initially, coarser afterwards) ----
dt_max = total_sim_time / n_coarse_steps;
dt_min = dt_max / 1000;
growth = (dt_max/dt_min)^(1/n_ramp_steps);

%% =========================================================================
%% SECTION 3: THERMAL MODEL OF MICRO-EDM   (paper Section 2.1, eqs. 1-5)
%% =========================================================================
% Empirical heat-flux-radius correlation (Yeo, Kurnia & Tan, 2007); current
% in Amps and pulse-on time in seconds, result in micrometers -- see note (2) above.
R_um = 2040 * I_curr^0.43 * Ton^0.44;        % [micrometers]
R    = R_um * 1e-6;                          % heat-flux radius [m]

q_heat = Cfrac*V*I_curr / (pi*R^2);          % heat flux [W/m^2], eq.(2)

tau = alpha_th*Ton / R^2;                    % non-dimensional time, eq.(5)

theta_peak = 2*sqrt(tau/pi)*(1 - exp(-1/(4*tau))) + erfc(1/(2*sqrt(tau))); % eq.(3)

T_peak = theta_peak * (Cfrac*V*I_curr) / (K_th*pi*R);  % dimensional peak T [K], eq.(4) rearranged

Ti = 0.5*(T_peak + Tb);   % initial (average) vapor temperature, per paper Sec. 2.1

fprintf('--- Thermal model results ---\n');
fprintf('Heat flux radius R      = %.3f um\n', R_um);
fprintf('Heat flux q             = %.3e W/m^2\n', q_heat);
fprintf('Peak temperature T_peak = %.1f K\n', T_peak);
fprintf('Initial vapor temp Ti   = %.1f K\n\n', Ti);

%% =========================================================================
%% SECTION 4-5: NUCLEATION-COAGULATION + COOLING-RATE MODEL
%%              (paper Sections 2.2-2.3, solved together per Fig.6)
%% =========================================================================

% ---- Initial conditions --------------------------------------------------
T  = Ti;
Ps_i = clausius_clapeyron(Ti, P0, Hlv, Rgas_p, Tb);
N1 = Ps_i / (kb*Ti);          % initial monomer number density [1/m^3]
Nk = zeros(n_nodes,1);        % number density in each particle-size node [1/m^3]

t  = 0;
dt = dt_min;

% ---- preallocate history --------------------------------------------------
maxsteps   = n_ramp_steps + n_coarse_steps + 10;
t_hist     = zeros(maxsteps,1);
Tt_hist    = zeros(maxsteps,1);
davg_hist  = zeros(maxsteps,1);
dstar_hist = NaN(maxsteps,1);
step = 0;

while t < total_sim_time && step < maxsteps
    step = step + 1;

    % ---- 4.1 Saturation state & critical diameter (eqs. 12-13) --------
    Ps = clausius_clapeyron(T, P0, Hlv, Rgas_p, Tb);
    Ns = Ps/(kb*T);
    S  = N1/Ns;

    if S > 1
        d_star = 4*sigma*v1 / (kb*T*log(S));       % eq.(12)
        v_star = (pi/6)*d_star^3;

        % ---- 4.2 Nucleation rate (eq. 14-15, classical nucleation theory) ----
        phi = sigma*a1/(kb*T);
        Jk  = Ns*N1*v1*sqrt(2*sigma/(pi*m1)) * exp( phi - 4*phi^3/(27*log(S)^2) ); % eq.(14)

        % locate bracketing node for v_star and the nucleation size operator eps_k, eq.(18-19)
        knuc = find(v_node >= v_star, 1, 'first');
        if isempty(knuc), knuc = n_nodes; end
        eps_k = min(v_star / v_node(knuc), 1);
    else
        d_star = NaN; v_star = 0; Jk = 0; knuc = []; eps_k = 0;
    end

    % ---- 4.3 Coagulation (eqs. 16, 20-21) ------------------------------
    beta  = const_beta * sqrt(6*kb*T/rho_p) * G;      % [n_nodes x n_nodes]
    Mmat  = beta .* (Nk*Nk');
    Birth = accumarray(idxK,   0.5*wK  .*Mmat(:), [n_nodes,1]) + ...
            accumarray(idxKp1, 0.5*wKp1.*Mmat(:), [n_nodes,1]);
    Death = Nk .* (beta*Nk);
    dNk_coag = Birth - Death;

    dNk_nuc = zeros(n_nodes,1);
    if ~isempty(knuc)
        dNk_nuc(knuc) = Jk*eps_k;
    end

    % ---- 4.4 Update particle-size distribution -------------------------
    Nk = max(Nk + dt*(dNk_nuc + dNk_coag), 0);
    N1 = max(N1 - dt*Jk*(v_star/v1), 0);

    % ---- 4.5 Average particle diameter (number-weighted volume mean) ---
    if sum(Nk) > 0
        v_avg = sum(Nk.*v_node)/sum(Nk);
        d_avg = (6*v_avg/pi)^(1/3);
    else
        d_avg = 0;
    end

    % ---- 4.6 Cooling-rate model (paper Section 2.2, eqs. 6-11) ---------
    if d_avg > 0
        d_p = d_avg;
    elseif ~isnan(d_star)
        d_p = d_star;
    else
        d_p = (6*v1/pi)^(1/3);
    end
    Pr    = mu_f*Cp_f/kf;
    Do    = kb*T / (3*pi*mu_f*d_p);           % eq.(8)
    v_vel = Do/lambda_mfp;                     % eq.(7)
    Re    = rho_f*v_vel*d_p/mu_f;
    Nu    = Re^2 * Pr^2;                       % eq.(6)
    h     = Nu*kf/d_p;                         % eq.(9)
    coolrate = 6*h/(rho_p*Cp_p*d_p) * (T - Tf);% eq.(11)
    T = max(T - coolrate*dt, Tf);

    % ---- store history --------------------------------------------------
    t_hist(step)     = t;
    Tt_hist(step)    = T;
    davg_hist(step)  = d_avg*1e9;      % nm
    dstar_hist(step) = d_star*1e9;     % nm (NaN once S<=1, i.e. nucleation has stopped)

    % ---- advance time & ramp time-step -----------------------------------
    t  = t + dt;
    dt = min(dt*growth, dt_max);
end

t_hist     = t_hist(1:step);
Tt_hist    = Tt_hist(1:step);
davg_hist  = davg_hist(1:step);
dstar_hist = dstar_hist(1:step);

fprintf('--- Simulation summary ---\n');
fprintf('Final temperature reached : %.1f K\n', Tt_hist(end));
fprintf('Final average diameter    : %.2f nm\n', davg_hist(end));
fprintf('Minimum critical diameter : %.3f nm\n\n', min(dstar_hist));

%% =========================================================================
%% SECTION 6: RESULTS & PLOTS
%% =========================================================================

figure('Name','Average particle diameter vs time');
plot(t_hist*1e6, davg_hist, 'b-', 'LineWidth', 1.8);
xlabel('Time (\mus)');
ylabel('Average particle diameter (nm)');
title(sprintf('Average nanoparticle diameter vs time  (V=%gV, I=%gA, T_{on}=%g\\mus)', ...
    V, I_curr, Ton*1e6));
grid on;

% NOTE: the three helper functions used above (clausius_clapeyron,
% build_coagulation_operator, build_collision_geometry) are defined in
% their own separate .m files in this same folder -- keep all four files
% together in one directory when you run this script.
