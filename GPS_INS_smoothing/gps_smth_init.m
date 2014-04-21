
function [nav] = gps_smth_init(data, p, start_time, yaw)
I3 = eye(3);
% --- initialize parameters ---

%p = gpsins_param();

% --- initialize measurement indices ---

nav.ii_imu = 1;

% --- initialize start time ---

nav.t_cpu = data.imu(nav.ii_imu).imu_tm;
nav.t_inertial = 0;

% --- initialize navigation state (needs improvement) ---

% GPS positioning result to provide initial position
rpos_ecef = LS_double_diff(data.gps(1), p.ecef_p_b, 1, 1);
%[rpos_ecef_2, delta_t_2] = LS_single_diff_L2(data.gps(1), p.ecef_p_b);
rpos_tp=p.R_ned2ecef' *( rpos_ecef - p.ecef_p_b);

roll = 0;
pitch = 0;
%yaw = 0;
%yaw = init_yaw(p, data.gps(2), data.gps(1));

xhat = gpsins_reset_xhat(nav.t_cpu, p.base_pos, roll, pitch, yaw); % change the start_time to nav.t_cpu
ic = randomcov(1000, I3*1^2, [0,0,0]); % generate 1000 random vector with covariance and mean [0,0,0]
ic(1,:) = rpos_tp;%[-18.1579647328271,5.3658874799492,-14.0794410885173]; 
xhat.r_tb_t = p.base_pos+ic(1,:)';%;%
xhat.v_tb_t = zeros(3,1);
xhat.ba_b = zeros(3,1);
xhat.bg_b = zeros(3,1);
%10*ic(1,:)'




% --- initialize error state (we can do far better here) ---

dx.Qd = [];
dx.dx_minus = zeros(p.X_STATES, 1);
dx.Pxx_minus = zeros(p.X_STATES, p.X_STATES);
dx.Pxx_minus(p.X_ANG, p.X_ANG) = I3*0.1^2;
dx.Pxx_minus(p.X_ANG(3), p.X_ANG(3)) = pi^2;
dx.Pxx_minus(p.X_VEL, p.X_VEL) = I3*0.5^2;
% use GPS positioning result to provide initial position
    dx.Pxx_minus(p.X_POS, p.X_POS) = diag([5,5,5].^2);
%dx.Pxx_minus(p.X_POS, p.X_POS) = I3*3^2;
dx.Pxx_minus(p.X_BA, p.X_BA) = I3*0.1^2;
dx.Pxx_minus(p.X_BG, p.X_BG) = I3*0.1^2;
dx.ydim = 0;

% --- create nav structure ---

nav.param = p;
nav.xhat = xhat;
nav.xhat0 = xhat;
nav.dx = dx;

