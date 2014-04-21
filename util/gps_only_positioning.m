clc;
clear;
close all;
%% options
Start_time = 17;
End_time = 157;
DGPS_ONLY = 1;
CODE_ONLY = 0;
%% constants
p = gpsins_param();

%% load data
load GPS_buff_log22.dat;
Buff = read_GPS_buff_log(GPS_buff_log22, Start_time, End_time);
buff_size = End_time-Start_time+1;
trajectory = zeros(buff_size,4);

for n = 1:buff_size; 
    meas_temp = Buff(n);
    rpos_tp=p.R_ned2ecef' *(LS_double_diff(meas_temp, p.ecef_p_b, DGPS_ONLY) - p.ecef_p_b);
    trajectory(n,1)  = meas_temp.imu_tm;
    trajectory(n,2:4)= rpos_tp; % save the position to trajectory      
end

fig = 2;
h = figure(fig);
fig = fig + 1;
subplot(3,1,1);
plot(trajectory(:,1), trajectory(:,2)-p.true_pos(1),'*');
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(trajectory(:,1), trajectory(:,3)-p.true_pos(2),'*');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(trajectory(:,1), trajectory(:,4)-p.true_pos(3),'*');
hold on;
grid on;
ylabel('z (m)');


