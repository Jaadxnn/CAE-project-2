function T_mean = mean_interface_temp(a2, L, Pz, Pxy)
%Inputs obtained from cutting simulations/experiments:
%a2-Deformed chip thickness mm
%L-Chip-tool contact length along the rake face m
%Pz-Main cutting force acting in the cutting direction N
%Pxy-Feed/thrust force acting normal to the cutting direction N


fprintf('=========================================================\n');
fprintf('  MEAN CHIP-TOOL INTERFACE TEMPERATURE\n');
fprintf('  Method: Singamneni (AUT Lecture Notes, 2024)\n');
fprintf('=========================================================\n\n');
%Process parameters
Vc = 200;
t = 2.5;
s0 = 0.3;

%Tool geometry
gamma0 = 0;

%Convert contact length to mm
L = L*1000;

%Material properties
K = 48;
rho = 8000;
C = 500;

%Assumed constants
Bi = 0.15;
theta_i = 30;

%Uncut chip geometry
phi_s = 90;
a1 = s0*sind(phi_s);
b1 = t/sind(phi_s);

%Chip compression ratio and shear angle
zeta = a2/a1;
beta = atan(cosd(gamma0)/(zeta-sind(gamma0)));
beta_deg = rad2deg(beta);

%Chip velocity
Vf = Vc/zeta;

%Cutting force components
F = abs(Pz)*sind(gamma0)+abs(Pxy)*cosd(gamma0);
N = abs(Pz)*cosd(gamma0)-abs(Pxy)*sind(gamma0);
mu = F / N;
Ps = abs(Pz)*cos(beta)-abs(Pxy)*sin(beta);

%Convert dimensions to metres
a2_m = a2*1e-3;
b1_m = b1*1e-3;

%Shear zone temperature
theta_s = (Ps*cosd(gamma0)) / ...
          (rho*C*b1_m*a2_m*cos(beta))* ...
          (1-Bi)+theta_i;

%Thermal parameter
L2 = (Vf*L*rho*C)/(4*K*1000*60);

%Chip-tool interface temperature rise
theta_t = (0.377*a2*F*Vf*1e3)/...
          (60*b1*K*sqrt(L2));

%Mean interface temperature
T_mean = theta_s + theta_t;

%Display results
fprintf('\nMean Interface Temperature Results\n');
fprintf('----------------------------------\n');
fprintf('zeta     = %.4f\n', zeta);
fprintf('beta     = %.2f deg\n', beta_deg);
fprintf('Vf       = %.4f m/min\n', Vf);
fprintf('F        = %.2f N\n', F);
fprintf('N        = %.2f N\n', N);
fprintf('mu       = %.4f\n', mu);
fprintf('Ps       = %.2f N\n', Ps);
fprintf('L2       = %.4f\n', L2);
fprintf('theta_s  = %.2f *C\n', theta_s);
fprintf('theta_t  = %.2f *C\n', theta_t);
fprintf('T_mean   = %.2f *C\n', T_mean);

end