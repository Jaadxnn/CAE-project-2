function T_mean = mean_interface_temp(a2, L, Pz, Pxy)
% =========================================================
%  MEAN CHIP-TOOL INTERFACE TEMPERATURE CALCULATOR
%  Method: Singamneni (AUT Lecture Notes, 2024)
%
%  T_mean = theta_s + theta_t
%
%  where:
%    theta_s = shear zone temperature rise
%    theta_t = mean interface temperature rise (friction)
%
%  Inputs (defined in USER INPUTS section below):
%    Vc    - Cutting speed        [m/min]
%    t     - Depth of cut         [mm]
%    s0    - Feed                 [mm/rev]
%    a2    - Chip thickness       [mm]
%    L     - Chip-tool contact    [mm]
%    Pz    - Main cutting force   [N]
%    Pxy   - Feed force           [N]
%    K     - Thermal conductivity [W/m.°C]
%    rho   - Density              [kg/m³]
%    C     - Specific heat        [J/kg.°C]
%    Bi    - Heat fraction into workpiece [-]
%    theta_i - Initial temperature [°C]
%    gamma0  - Rake angle          [deg]
% =========================================================


fprintf('=========================================================\n');
fprintf('  MEAN CHIP-TOOL INTERFACE TEMPERATURE\n');
fprintf('  Method: Singamneni (AUT Lecture Notes, 2024)\n');
fprintf('=========================================================\n\n');

% ---------------------------------------------------------
%  USER INPUTS
% ---------------------------------------------------------

% --- Process Parameters ---
Vc      = 200;      % Cutting speed        [m/min]
t       = 2.5;      % Depth of cut         [mm]
s0      = 0.3;      % Feed                 [mm/rev]

% --- Cutting Tool Geometry ---
gamma0  = 0;        % Rake angle           [deg]

% --- Measured / Given Responses ---
%{
a2      = 0.8;      % Chip thickness       [mm]   (bevelled tool)
L       = 0.30;     % Chip-tool contact length [mm] (diagonal for bevelled)
Pz      = 1050;     % Main cutting force   [N]
Pxy     = 750;      % Feed force           [N]
%}

L = L*1000;
% --- Work Material Properties (Mild Steel) ---
K       = 48;       % Thermal conductivity [W/m.°C]
rho     = 8000;     % Density              [kg/m³]
C       = 500;      % Specific heat        [J/kg.°C]

% --- Assumed Constants ---
Bi      = 0.15;     % Fraction of shear heat into workpiece [-]
theta_i = 30;       % Initial (ambient) temperature        [°C]

% ---------------------------------------------------------
%  STEP 1: UNCUT CHIP THICKNESS & WIDTH OF CUT
%  a1 = s0 * sin(phi_s)  where phi_s = 90 deg (orthogonal)
%  b1 = t / sin(phi_s)
% ---------------------------------------------------------
phi_s = 90;                         % cutting edge angle [deg] (orthogonal)
a1    = s0 * sind(phi_s);           % uncut chip thickness [mm]
b1    = t  / sind(phi_s);           % width of cut [mm]

fprintf('--- Step 1: Geometry ---\n');
fprintf('  Uncut chip thickness  a1 = s0 * sin(phi_s) = %.4f mm\n', a1);
fprintf('  Width of cut          b1 = t / sin(phi_s)  = %.4f mm\n\n', b1);

% ---------------------------------------------------------
%  STEP 2: CHIP COMPRESSION RATIO & SHEAR ANGLE
% ---------------------------------------------------------
zeta   = a2 / a1;                               % chip reduction coefficient
gamma0_rad = deg2rad(gamma0);
beta   = atan(cosd(gamma0) / (zeta - sind(gamma0)));  % shear angle [rad]
beta_deg = rad2deg(beta);

fprintf('--- Step 2: Kinematics ---\n');
fprintf('  Chip reduction coeff  zeta = a2/a1 = %.4f\n', zeta);
fprintf('  Shear angle           beta = %.4f deg\n\n', beta_deg);

% ---------------------------------------------------------
%  STEP 3: CHIP VELOCITY
% ---------------------------------------------------------
Vf = Vc / zeta;     % chip velocity [m/min]

fprintf('--- Step 3: Chip Velocity ---\n');
fprintf('  Chip velocity  Vf = Vc / zeta = %.4f m/min\n\n', Vf);

% ---------------------------------------------------------
%  STEP 4: FORCE COMPONENTS (Merchant's Circle, rake = 0)
%  F  = Pz*sin(gamma0) + Pxy*cos(gamma0)
%  N  = Pz*cos(gamma0) - Pxy*sin(gamma0)
%  Ps = Pz*cos(beta)   - Pxy*sin(beta)
% ---------------------------------------------------------
F  = Pz * sind(gamma0) + abs(Pxy) * cosd(gamma0);   % friction force  [N]
N  = Pz * cosd(gamma0) - abs(Pxy) * sind(gamma0);   % normal force    [N]
mu = F / N;                                     % friction coefficient
Ps = Pz * cos(beta) - abs(Pxy) * sin(beta);          % shear force     [N]

fprintf('--- Step 4: Forces ---\n');
fprintf('  Friction force  F  = %.2f N\n',   F);
fprintf('  Normal force    N  = %.2f N\n',   N);
fprintf('  Friction coeff  mu = %.4f\n',     mu);
fprintf('  Shear force     Ps = %.2f N\n\n', Ps);

% ---------------------------------------------------------
%  STEP 5: SHEAR ZONE TEMPERATURE RISE (theta_s)
%
%  theta_s = [ Ps * cos(gamma0) / (rho*C * b1[m] * a2[m] * cos(beta)) ]
%            * (1 - Bi) + theta_i
% ---------------------------------------------------------
a2_m = a2 * 1e-3;      % convert mm -> m
b1_m = b1 * 1e-3;      % convert mm -> m

theta_s = (Ps * cosd(gamma0)) / ...
          (rho * C * b1_m * a2_m * cos(beta)) ...
          * (1 - Bi) + theta_i;

fprintf('--- Step 5: Shear Zone Temperature ---\n');
fprintf('  theta_s = [Ps*cos(gamma0) / (rho*C*b1*a2*cos(beta))] * (1-Bi) + theta_i\n');
fprintf('  theta_s = %.2f degC\n\n', theta_s);

% ---------------------------------------------------------
%  STEP 6: THERMAL PARAMETER L2
%
%  L2 = Vf[m/min] * L[mm] * rho * C / (4 * K * 1000 * 60)
% ---------------------------------------------------------
L2 = (Vf * L * rho * C) / (4 * K * 1000 * 60);

fprintf('--- Step 6: Thermal Parameter ---\n');
fprintf('  L2 = Vf*L*rho*C / (4*K*1000*60)\n');
fprintf('  L2 = %.4f\n\n', L2);

% ---------------------------------------------------------
%  STEP 7: INTERFACE TEMPERATURE RISE (theta_t)
%
%  theta_t = 0.377 * a2[mm] * F * Vf[m/min] * 1e3
%            / (60 * b1[mm] * K * sqrt(L2))
% ---------------------------------------------------------
theta_t = (0.377 * a2 * F * Vf * 1e3) / ...
          (60 * b1 * K * sqrt(L2));

fprintf('--- Step 7: Interface Temperature Rise ---\n');
fprintf('  theta_t = 0.377 * a2 * F * Vf * 1e3 / (60 * b1 * K * sqrt(L2))\n');
fprintf('  theta_t = %.2f degC\n\n', theta_t);

% ---------------------------------------------------------
%  STEP 8: MEAN CHIP-TOOL INTERFACE TEMPERATURE
% ---------------------------------------------------------
T_mean = theta_s + theta_t;

fprintf('=========================================================\n');
fprintf('  RESULTS SUMMARY\n');
fprintf('=========================================================\n');
fprintf('  Cutting speed          Vc    = %.1f m/min\n',  Vc);
fprintf('  Feed                   s0    = %.2f mm/rev\n', s0);
fprintf('  Depth of cut           t     = %.1f mm\n',    t);
fprintf('  Chip thickness         a2    = %.2f mm\n',    a2);
fprintf('  Contact length         L     = %.2f mm\n',    L);
fprintf('  Main cutting force     Pz    = %.0f N\n',     Pz);
fprintf('  Feed force             Pxy   = %.0f N\n',     Pxy);
fprintf('---------------------------------------------------------\n');
fprintf('  Chip reduction coeff   zeta  = %.4f\n',       zeta);
fprintf('  Shear angle            beta  = %.2f deg\n',   beta_deg);
fprintf('  Chip velocity          Vf    = %.4f m/min\n', Vf);
fprintf('  Thermal parameter      L2    = %.4f\n',       L2);
fprintf('---------------------------------------------------------\n');
fprintf('  Shear zone temp        theta_s = %7.2f degC\n', theta_s);
fprintf('  Interface temp rise    theta_t = %7.2f degC\n', theta_t);
fprintf('---------------------------------------------------------\n');
fprintf('  MEAN INTERFACE TEMP    T_mean  = %7.2f degC\n', T_mean);
fprintf('=========================================================\n');

end