clc;
%clear;
%close all;

%% Add routine path
addpath('../');

%% options
Start_time = 10;
End_time = 330;
DGPS_ONLY = 0;
CODE_ONLY = 0;
%% constants
p = gpsins_param();

%% load data
% gps_data = load('../Data_set/GPS_buff_log22.dat');
Buff = read_GPS_buff_log(gps_data, Start_time, End_time);
buff_size = End_time-Start_time+1;
trajectory = zeros(buff_size,4);

for n = 1:buff_size; 
    meas_temp = Buff(n);
    rpos_tp=p.R_ned2ecef' *(LS_double_diff(meas_temp, p.ecef_p_b, DGPS_ONLY) - p.ecef_p_b);
    trajectory(n,1)  = meas_temp.imu_tm;
    trajectory(n,2:4)= rpos_tp; % save the position to trajectory      
end

fig = 1;
h = figure(fig);
fig = fig + 1;
subplot(3,1,1);
plot(trajectory(:,1), trajectory(:,2),'--+');
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(trajectory(:,1), trajectory(:,3),'--+');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(trajectory(:,1), trajectory(:,4),'--+');
hold on;
grid on;
ylabel('z (m)');

fig = 101;
h = figure(fig);
plane_fig = fig;
fig = fig + 1;
plot(trajectory(:,3), trajectory(:,2),'--+');%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '-g'
hold on; grid on;
title(['Plane trajectory'] );
ylabel('North, meter');
xlabel('East, meter');
