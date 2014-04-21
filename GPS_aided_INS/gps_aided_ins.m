% Author: Yiming Chen, yimingchen86@gmail.com
% This is a script for reprocessing the GPS and IMU data to implement
% GPS aided INS navigation system with EKF. Two test data sets are included:
% one is xxx_buff_log_veh.dat logged on moving vehicle and another is
% xxx_buff_log_veh.dat logged from the fixed antenna sitting still.

clc;
%clear;
%close all;

%% Add routine path
addpath('../');
addpath('../util');

%% Options
Start_time = 60;
End_time = 80;
DGPS = 1;  % only use single differenced L1 code
LOOSE = 1; % for loosely coupling

%% Load data
%  imu_data = load('../Data_set/IMU_buff_log_veh.dat');
%  gps_data = load('../Data_set/GPS_buff_log_veh.dat');
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

%% Process algorithm
% initialize navigation state
nav = gpsins_init(data, p, Start_time);

for ii_imu = 1:N
	% perform one navigation iteration
	[nav, out] = gpsins_step(p, nav, ii_imu, data, out, DGPS, LOOSE);
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

fig = 100;
h = figure(fig);
fig = fig + 1;
subplot(3,1,1);
plot(out.t, out.r_tb_t(:,1));%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '-g'
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2));%, out.t,[3*sqrt(out.Pxx(:,8)),-3*sqrt(out.Pxx(:,8))], '-g'
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3));%, out.t,[3*sqrt(out.Pxx(:,9)),-3*sqrt(out.Pxx(:,9))], '-g'
hold on;
grid on;
ylabel('z (m)');

h = figure(fig);
plane_fig = fig;
fig = fig + 1;
plot(out.r_tb_t(:,2), out.r_tb_t(:,1));%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '-g'
hold on; grid on;
title(['Plane trajectory'] );
ylabel('North, meter');
xlabel('East, meter');


