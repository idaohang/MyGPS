% This script is for gps-ins smoothing only. No near-real-time approach
function [N_error1,N_error2, N_error3, num_sats, X_hat1, X_hat2, X_hat3, out1, out2, out3, out0, mils_res] = gps_smth_ins(Start_time,Buff_Size,state_out)
%clc;
%clear;
%close all;

global imu_data;
global gps_data;
r2d = 180/pi;

%% Options
% Start_time = 105;
% Buff_Size = 10;
End_time = Start_time+Buff_Size;
global DGPS;
global Warm_Up;
global warm_time;
global Dual_Freq;
global Stationary;
LOOSE = 0;     % for loosely coupling
Warm_Up = 0;
warm_time = 0;

%% Load data
% imu_data = load('../Data_set/IMU_buff_log_veh.dat');
% gps_data = load('../Data_set/GPS_buff_log_veh.dat');
data.imu = read_IMU_buff_log(imu_data, Start_time, End_time);
[~,m] = size(gps_data);
if m>1525
    data.gps = read_GPS_buff_log_new(gps_data, Start_time, End_time);
else
    data.gps = read_GPS_buff_log(gps_data, Start_time, End_time);
end
%state_out = read_smooth_result(state_out, Start_time, End_time);

% load('../Data_set/python_smooth.mat');
% % trim data for plotting
% i = 1;
% while ( state_out(i,1) < Start_time )
%     i=i+1;
% end
% j=i+1;
% while ( state_out(j,1)<Start_time+ Buff_Size-0.5)
%     j=j+1;
% end
% state_out = state_out(i:j,:);

%% Inject Outlier
%data = eraim_outlier_injection(data);

%% Simulation setup
% simulation parameters
p = gpsins_param();

% number of IMU samples to process
N = size(data.imu,1);

% allocate memory for data analysis
out.t = ones(N,1)*NaN;
out.r_tb_t = ones(N,3)*NaN;
out.v_tb_t = ones(N,3)*NaN;
out.v_tb_b = ones(N,3)*NaN;
out.ba_b = ones(N,3)*NaN;
out.bg_b = ones(N,3)*NaN;
out.alphahat = ones(N,3)*NaN;
out.Pxx = ones(N,p.X_STATES)*NaN;

% process algorithm
if isempty(state_out) % for stationary case
    yaw = 0;
else
    yaw = state_out(1,7);
end
nav = gps_smth_init(data, p, Start_time, yaw);
nav = gps_smth_buff_init(nav, Buff_Size);

N_error1 = [];
N_error2 = []; 
num_sats = [];
X_hat1   = [];
X_hat2   = []; 
out1     = [];
out2     = [];
out3     = [];

for ii_imu = 1:N
    % perform one navigation iteration
    [nav, out] = gps_smth_step(p, nav, ii_imu, data, out, LOOSE);
    
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
    
    if (nav.buff.kdim==Buff_Size)
        break
    end
end

% figure(3);
% plot(out.t, limit_pi(out.alphahat)*r2d, '.r');
% hold on;
% grid on;

prn_list = merge_prn_list(data);
num_sats = length(prn_list);

if isempty(prn_list)
    [N_error1] = NaN;
    [N_error2] = NaN;
    display('Empty Prn list!!! No valid sateliites!!!')
    return;
end

Xhat_0 = [nav.buff.xhat_k0;nav.buff.xhat_k];
[smth_out, X_hat] = gps_smoothing(nav,data,prn_list);

Phase_meas = validate_phase_meas(nav, data, X_hat);
if length(Phase_meas.Nhat) > 6
    %[smth_out, X_hat3, N_hat3] = gps_smoothing_int(nav,data,X_hat,Phase_meas);
    X_hat3 = X_hat; % FIXME
    N_hat3 = Phase_meas.Nhat;
    [smth_out, X_hat1, N_hat1] = gps_smoothing_shift(nav, data, X_hat3, Phase_meas);
    %Phase_meas.Nhat = N_hat1;
    %[rt, Jt, gt, mils_res] = gps_smth_res_J_mils(nav, data, X_hat1, Phase_meas);
    %X_hat2 = X_hat1; % FIXME
    %N_hat2 = N_hat1;
    [smth_out, X_hat2, N_hat2, mils_res] = gps_smoothing_miles(nav,data,X_hat3,Phase_meas);
else
    X_hat1 = X_hat; % FIXME
    N_hat1 = Phase_meas.Nhat;
    X_hat2 = X_hat; % FIXME
    N_hat2 = Phase_meas.Nhat;
    X_hat3 = X_hat; % FIXME
    N_hat3 = Phase_meas.Nhat;
    mils_res = [];
end

% Propagate whole trajectory
out0 = prop_whole_traj(nav, data, X_hat);
out0 = [out0.t, out0.r_tb_t];
out1 = prop_whole_traj(nav, data, X_hat1);
out1= [out1.t, out1.r_tb_t];
out2 = prop_whole_traj(nav, data, X_hat2);
out2 = [out2.t, out2.r_tb_t];
out3 = prop_whole_traj(nav, data, X_hat3);
out3 = [out3.t, out3.r_tb_t];

% m = size(state_out);
% out11 = [out1.t([1:2:2*m],:), out1.r_tb_t([1:2:2*m],:)];
% error1 = mean(state_out(:,2:4) - out11(:,2:4));
% var1 = var(state_out(:,2:4) - out11(:,2:4));
% out12 = [out2.t([1:2:2*m],:), out2.r_tb_t([1:2:2*m],:)];
% error2 = mean(state_out(:,2:4) - out12(:,2:4));
% var2 = var(state_out(:,2:4) - out12(:,2:4));
% display(['Positioning error1 mean is ', num2str(error1)] );
% display(['Positioning error2 mean ', num2str(error2)] );
% display(['Positioning error1 variance is ', num2str(var1)] );
% display(['Positioning error2 variance is ', num2str(var2)] );

N_hat = check_integer(Phase_meas, data.gps(end-1,1));
N_error1 =  N_hat + N_hat1 ;
display(['Integer error norm1 is ', num2str( sqrt(norm(N_error1)^2/length(N_error1)) )] );
N_error2 =  N_hat + N_hat2;
display(['Integer error norm2 is ', num2str( sqrt(norm(N_error2)^2/length(N_error2)) )] );
N_error3 =  N_hat + N_hat3;
display(['Integer error norm3 is ', num2str( sqrt(norm(N_error3)^2/length(N_error3)) )] );
%plot_d_Xhat(X_hat, X_hat1);O

% fig = 1;
% h = figure(fig);
% fig = fig + 1;
% subplot(3,1,1);
% plot(out.t, out.r_tb_t(:,1)-p.base_pos(1), '.k');%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '--g'
% title(['NED position'] );
% hold on;
% grid on;
% ylabel('x (m)');
% subplot(3,1,2);
% plot(out.t, out.r_tb_t(:,2)-p.base_pos(2), '.k');%, out.t,[3*sqrt(out.Pxx(:,8)),-3*sqrt(out.Pxx(:,8))], '--g'
% hold on;
% grid on;
% ylabel('y (m)');
% subplot(3,1,3);
% plot(out.t, out.r_tb_t(:,3)-p.base_pos(3), '.k');%, out.t,[3*sqrt(out.Pxx(:,9)),-3*sqrt(out.Pxx(:,9))], '--g'
% hold on;
% grid on;
% ylabel('z (m)');

%Plot Smoothed result from python code
% if ~isempty(state_out)
%     figure(1);
%     subplot(3,1,1);
%     plot(state_out(:,1), state_out(:,2), '.k', out2(:,1), out2(:,2), '.g' );
%     title(['NED position'] );
%     hold on;
%     grid on;
%     ylabel('x (m)');
%     subplot(3,1,2);
%     plot(state_out(:,1), state_out(:,3), '.k', out2(:,1), out2(:,3), '.g' );
%     hold on;
%     grid on;
%     ylabel('y (m)');
%     subplot(3,1,3);
%     plot(state_out(:,1), state_out(:,4), '.k', out2(:,1), out2(:,4), '.g' );
%     hold on;
%     grid on;
%     ylabel('z (m)');
% 
%     legend('smoothing 1', 'smoothing 2', 'RTK EKF')
% 
%     % figure(2);
%     % subplot(3,1,1);
%     % plot(pos_ins(:,1), pos_ins(:,5), '.k');
%     % hold on;
%     % grid on;
%     % subplot(3,1,2);
%     % plot(pos_ins(:,1), pos_ins(:,6), '.k');
%     % hold on;
%     % grid on;
%     % subplot(3,1,3);
%     % plot(pos_ins(:,1), pos_ins(:,7), '.k');
%     % hold on;
%     % grid on;
% 
%     figure(3);
%     plot(state_out(:,1), state_out(:,7), '.k');
%     hold on;
%     grid on;
% 
%     figure(4);
%     plot(state_out(:,3), state_out(:,2), '.k');
%     legend('Tangent Trajectory')
%     hold on;
%     grid on;
% end

%% Plot RTK EKF result from C code
% figure(1);
% subplot(3,1,1);
% plot(pos_ins(:,1), pos_ins(:,2), '.k');
% title(['NED position'] );
% hold on;
% grid on;
% ylabel('x (m)');
% subplot(3,1,2);
% plot(pos_ins(:,1), pos_ins(:,3), '.k');
% hold on;
% grid on;
% ylabel('y (m)');
% subplot(3,1,3);
% plot(pos_ins(:,1), pos_ins(:,4), '.k');
% hold on;
% grid on;
% ylabel('z (m)');
% 
% legend('smoothing 1', 'smoothing 2', 'RTK EKF')
% 
% figure(2);
% subplot(3,1,1);
% plot(pos_ins(:,1), pos_ins(:,5), '.k');
% hold on;
% grid on;
% subplot(3,1,2);
% plot(pos_ins(:,1), pos_ins(:,6), '.k');
% hold on;
% grid on;
% subplot(3,1,3);
% plot(pos_ins(:,1), pos_ins(:,7), '.k');
% hold on;
% grid on;
% 
% figure(3);
% plot(pos_ins(:,1), pos_ins(:,10), '.k');
% hold on;
% grid on;
% 
% figure(4);
% plot(pos_ins(:,3), pos_ins(:,2), '.k');
% hold on;
% grid on;

%% Plot residuals
% h = figure(fig);
% fig = fig + 1;
% % for k = 1:size(smth_out.res.t)
% %     plot(smth_out.res.t(k), smth_out.res.r{k}, '+');%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '--g'
% %     hold on;
% % end
% % Res = cell2mat(smth_out.res.r);
% plot(smth_out.res.t, cell2mat(smth_out.res.r), '+--');%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '--g'
% hold on;
% grid on;
% title(['Residual'] );
% ylabel('residual (m)');


% h = figure(fig);
% fig = fig + 1;
% %plot(out.residual.t, out.residual.dy, '+--');%, out.t,[3*sqrt(out.Pxx(:,7)),-3*sqrt(out.Pxx(:,7))], '--g'
% hist(out.residual.dy(10:end,1));
% hold on;
% grid on;
% title(['Residual'] );
% ylabel('residual (m)');

% r2d = 180/pi;
% h = figure(3);
% color = get(gca,'ColorOrder');
% plot(out.t, limit_pi(out.alphahat)*r2d,'--*b');
% title(['Yaw angle'] );
% hold on;
% grid on;
% ylabel('yaw (deg)');


