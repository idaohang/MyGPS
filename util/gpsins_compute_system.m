
function [Phi, Qd, Q, G1] = gpsins_compute_system(p, xhat, ya_b, yg_b, dt)

[F, G, Q] = compute_system(p, xhat, ya_b, yg_b);
[Phi, G1] = compute_discrete(F, G, dt);
% Qd = G*Q*G';  % DZ: Qd computed by this equation would always be singular
                %    since G is a tall matrix. This approximation assumes Phi
                %    = I5, which is very rough.
Qd = G*Q*G'*dt + 1/2*F*G*Q*G'*(dt)^2 + 1/2*G*Q*G'*F'*(dt)^2 + 1/3*F*G*Q*G'*F'*(dt)^3;
% DZ: Using the above equation to compute Qd instead, where we assume Phi = I +
%    F*dt. This approximation would be more accurate and guarantee Qd is
%    nonsingular.


%compute continous time system
function [F, G, Q] = compute_system(p, xhat, ya_b, yg_b)

% allocate memory
F = zeros(p.X_STATES);
G = zeros(p.X_STATES, p.U_INPUTS);
Q = zeros(p.U_INPUTS);
I3 = eye(3);

% compute inertial quantities
omega_ib_b = yg_b - xhat.bg_b;
omega_tb_b = omega_ib_b - xhat.R_b2t'*p.omega_it_t;
f_ib_b = ya_b - xhat.ba_b;

% attitude error state
F(p.X_ANG, p.X_BG) = -I3;
F(p.X_ANG, p.X_ANG) = -vcross(omega_tb_b);
G(p.X_ANG, p.U_NG) = -I3;

% velocity error state
F(p.X_VEL, p.X_BA) = -xhat.R_b2t;
F(p.X_VEL, p.X_ANG) = -xhat.R_b2t*vcross(f_ib_b);
F(p.X_VEL, p.X_VEL) = -2*vcross(p.omega_it_t);
G(p.X_VEL, p.U_NA) = -xhat.R_b2t;

% position error state
F(p.X_POS, p.X_VEL) = I3;

% accelerometer bias error state
G(p.X_BA, p.U_NBA) = I3;

% gyro bias error state
G(p.X_BG, p.U_NBG) = I3;

% noise
Q(p.U_NA, p.U_NA) = I3*p.sigma_na*p.sigma_na';
Q(p.U_NBA, p.U_NBA) = I3*p.sigma_nba*p.sigma_nba';
Q(p.U_NG, p.U_NG) = I3*p.sigma_ng*p.sigma_ng';
Q(p.U_NBG, p.U_NBG) = I3*p.sigma_nbg*p.sigma_nbg';

% compute discrete time system
function [Phi, G] = compute_discrete(F, G, dt)
Phi = eye(size(F)) + F*dt + (F*dt)^2/2;
G = G*sqrt(dt); % Qd = G*Q*G'

%% Old
% function [F, G, Q] = compute_system(p, xhat, ya_b, yg_b)
% 
% % allocate memory
% F = zeros(p.X_STATES);
% G = zeros(p.X_STATES, p.U_INPUTS);
% Q = zeros(p.U_INPUTS);
% I3 = eye(3);
% 
% % compute inertial quantities
% omega_ib_b = yg_b - xhat.bg_b;
% omega_tb_b = omega_ib_b - xhat.R_b2t'*p.omega_it_t;
% f_ib_b = ya_b - xhat.ba_b;
% 
% % attitude error state
% F(p.X_ANG, p.X_BG) = -I3;
% F(p.X_ANG, p.X_ANG) = -vcross(omega_tb_b);
% G(p.X_ANG, p.U_NG) = -I3;
% 
% % velocity error state
% F(p.X_VEL, p.X_BA) = -xhat.R_b2t;
% F(p.X_VEL, p.X_ANG) = -xhat.R_b2t*vcross(f_ib_b);
% F(p.X_VEL, p.X_VEL) = -2*vcross(p.omega_it_t);
% G(p.X_VEL, p.U_NA) = -xhat.R_b2t;
% 
% % position error state
% F(p.X_POS, p.X_VEL) = I3;
% 
% % accelerometer bias error state
% G(p.X_BA, p.U_NBA) = I3;
% 
% % gyro bias error state
% G(p.X_BG, p.U_NBG) = I3;
% 
% % noise
% Q(p.U_NA, p.U_NA) = I3*p.sigma_na*p.sigma_na';
% Q(p.U_NBA, p.U_NBA) = I3*p.sigma_nba*p.sigma_nba';
% Q(p.U_NG, p.U_NG) = I3*p.sigma_ng*p.sigma_ng';
% Q(p.U_NBG, p.U_NBG) = I3*p.sigma_nbg*p.sigma_nbg';
% 
% 
% % compute discrete time system
% function [Phi, G] = compute_discrete(F, G, dt)
% Phi = eye(size(F)) + F*dt + (F*dt)^2/2;
% G = G*sqrt(dt); % Qd = G*Q*G'

