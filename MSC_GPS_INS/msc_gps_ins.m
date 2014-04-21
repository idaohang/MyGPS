clc;
clear;
%close all;

%% Add routine path
addpath('../');

%% Options
Start_time = 25;%17;
End_time = 35;%57;
K = 10; % buffer length 

%% Load data
imu_data = load('../Data_set/IMU_buff_log22.dat');
gps_data = load('../Data_set/GPS_buff_log22.dat');
data.imu = read_IMU_buff_log(imu_data, Start_time, End_time);
data.gps = read_GPS_buff_log(gps_data, Start_time, End_time);

%% Simulation setup
% simulation parameters
p = gpsins_param();
% number of IMU samples to process
N = size(data.imu,1);
% allocate memory for data analysis
out.t = zeros(N,1);
out.r_tb_t = zeros(N,3);
out.v_tb_t = zeros(N,3);
out.v_tb_b = zeros(N,3);
out.ba_b = zeros(N,3);
out.bg_b = zeros(N,3);
out.alphahat = zeros(N,1);
out.Pxx = zeros(N,p.X_STATES);

out.t_msc = [];
out.r_tb_t_msc = [];

%% Process algorithm
% initialize navigation state
nav = msc_gps_init(data, p, Start_time);
buff_size = End_time-Start_time;
nav = msc_gps_buff_init(nav, K);

count = 0;
for ii_imu = 1:N
	% perform one navigation iteration
	[nav, out, count] = msc_gps_step(p, nav, ii_imu, data, out, K, count);
	% log data
	out.t(ii_imu,:) = nav.xhat.t;
	out.r_tb_t(ii_imu,:) = (nav.xhat.r_tb_t)';
	out.v_tb_t(ii_imu,:) = (nav.xhat.v_tb_t)';
	out.v_tb_b(ii_imu,:) = (nav.xhat.R_b2t'*nav.xhat.v_tb_t)';
    out.R_b2t(:,:,ii_imu)  = nav.xhat.R_b2t;
	out.ba_b(ii_imu,:) = nav.xhat.ba_b';
	out.bg_b(ii_imu,:) = nav.xhat.bg_b';
	out.Pxx(ii_imu,:) = diag(nav.dx.Pxx_minus);
	xt = nav.xhat.R_b2t*[1;0;0];
	out.alphahat(ii_imu,:) = atan2(xt(2),xt(1));

	if find(isnan(nav.dx.Pxx_minus))
		break
	end
end

fig = 2;
h = figure(fig);
fig = fig + 1;
subplot(3,1,1);
plot(out.t, out.r_tb_t(:,1)-p.base_pos(1), '-r', out.t_msc, out.r_tb_t_msc(:,1)-p.base_pos(1), '-go', out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '-b');

title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2)-p.base_pos(2), '-r', out.t_msc, out.r_tb_t_msc(:,2)-p.base_pos(2), '-go', out.t,[3*sqrt(out.Pxx(:,8)),-3*sqrt(out.Pxx(:,8))], '-b');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3)-p.base_pos(3), '-r', out.t_msc, out.r_tb_t_msc(:,3)-p.base_pos(3), '-go', out.t,[3*sqrt(out.Pxx(:,9)),-3*sqrt(out.Pxx(:,9))], '-b');
hold on;
grid on;
ylabel('z (m)');


