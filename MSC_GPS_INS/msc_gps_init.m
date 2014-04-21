
function [nav] = msc_gps_init(data, p, start_time)
I3 = eye(3);
% --- initialize parameters ---

%p = gpsins_param();

% --- initialize measurement indices ---

nav.ii_imu = 1;

% --- initialize start time ---

nav.t_cpu = data.imu(nav.ii_imu).imu_tm;
nav.t_inertial = 0;

% --- initialize navigation state (needs improvement) ---

xhat = gpsins_reset_xhat(start_time, p.base_pos);
% Later, use GPS positioning result to provide initial position 
 ic = randomcov(1000, I3*1^2, [0,0,0]); % generate 1000 random vector with covariance and mean [0,0,0]
ic(1,:) = [0,0,0]; %[2823.2,-908.25,57.44];+rand(1,3)
xhat.r_tb_t = p.base_pos+ic(1,:)';%;%zeros(3,1);
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
nav.dx = dx;

