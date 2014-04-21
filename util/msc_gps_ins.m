clc;
clear;
close all;
%% Options
Start_time = 17;
End_time = 57;
K = 30; % buffer length 

%% Load data
load ./Data_set/IMU_buff_log22.dat;
load ./Data_set/GPS_buff_log22.dat;
data.imu = read_IMU_buff_log(IMU_buff_log22, Start_time, End_time);
data.gps = read_GPS_buff_log(GPS_buff_log22, Start_time, End_time);

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

%% Process algorithm
% initialize navigation state
nav = gpsins_init(data, p, Start_time);
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
plot(out.t, out.r_tb_t(:,1)-p.true_pos(1), out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '-g');
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2)-p.true_pos(2), out.t,[3*sqrt(out.Pxx(:,8)),-3*sqrt(out.Pxx(:,8))], '-g');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3)-p.true_pos(3), out.t,[3*sqrt(out.Pxx(:,9)),-3*sqrt(out.Pxx(:,9))], '-g');
hold on;
grid on;
ylabel('z (m)');


