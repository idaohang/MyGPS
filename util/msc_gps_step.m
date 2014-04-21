
function [nav, out, count] = msc_gps_step(p, nav, ii_imu, data, out, K, count)

nav.dx.applied = 0;
nav.dx.ydim = 0;
N = nav.param.X_STATES;

% sample inertial sensor
nav.ii_imu = ii_imu;
nav.t_cpu = data.imu(ii_imu).imu_tm;
if ii_imu > 1
    t_cpu_1 = data.imu(ii_imu-1).imu_tm;
else
    t_cpu_1 = data.imu(ii_imu).imu_tm - data.imu(ii_imu).dt; %0.0050;
end
dt = data.imu(ii_imu).dt;
ya_b = data.imu(max(ii_imu-1,1)).accel_meas;
yg_b = data.imu(max(ii_imu-1,1)).gyro_meas;
%  ya_b = data.imu(ii_imu).accel_meas;
%  yg_b = data.imu(ii_imu).gyro_meas;

% propagate navigation state
[nav.xhat] = gpsins_propagate_xhat(nav.param, nav.xhat, ya_b, yg_b, dt);

% propagate error state about trajectory
[nav.dx, Phi_i, Qd_i] = gpsins_propagate_dx(nav.param, nav.dx, nav.xhat, ya_b, yg_b, dt);

%sample GPS measurement
ii_gps = gpsins_measurement_indices(data.gps, nav.t_cpu, t_cpu_1);
if ii_gps>0
    %     if count<K % if buffer not full
    %         pos_ecef = p.ecef_p_b + p.R_ned2ecef*nav.xhat.r_tb_t;
    %         [H, R, dy] = gpsins_gps_output(nav.param, pos_ecef, data.gps(ii_gps), 1); % always use DGPS not RTK
    %         [nav.dx, out] = gpsins_correct_dx(nav.dx, H, R, dy, out); % GPS aiding
    %         %count=count+1; % increment the buffer count
    %
    %         % correct navigation state
    %         if norm(nav.dx.dx_minus) > 0
    %             [nav.xhat] = gpsins_correct_xhat(nav.param, nav.xhat, nav.dx.dx_minus);
    %             if nav.dx.applied
    %                 count=count+1; % increment the buffer count
    %                 nav.Xhat(count) = nav.xhat;
    %                 nav.dx_buff(count) = nav.dx;
    %                 nav.Pxx = blkdiag(nav.Pxx,nav.dx.Pxx_minus);
    %             end
    %             nav.dx.dx_minus(:) = 0;
    %         end
    %
    %     else % count = K, buffer full
    %         pos_ecef = p.ecef_p_b + p.R_ned2ecef*nav.xhat.r_tb_t;
    %         [H, R, dy] = gpsins_gps_output(nav.param, pos_ecef, data.gps(ii_gps), 1); % always use DGPS not RTK
    %         [nav.dx, out] = gpsins_correct_dx(nav.dx, H, R, dy, out); % GPS aiding
    %
    %         nav.Xhat = [nav.Xhat(2:end);nav.xhat]; % pop xhat buffer
    %         Pos_ecef = zeros(K*3,1);
    %         for n=1:K
    %             Pos_ecef(3*n-2:3*n) = p.ecef_p_b + p.R_ned2ecef*nav.Xhat(n).r_tb_t;
    %         end
    %         nav.dx_buff = [nav.dx_buff(2:end), nav.dx]; % pop dx buffer
    %         nav.Pxx = blkdiag(nav.Pxx(N+1:end,N+1:end),nav.dx.Pxx_minus);
    %
    %         gps_plot_Xhat(p, nav.Xhat, nav.dx_buff, K, 100);
    %
    %         [H, R, dy] = msc_gps_ins_output(nav.param, Pos_ecef, data.gps(ii_gps-K+1:ii_gps), K);
    %         [nav] = msc_gps_correct_Dx(nav, H, R, dy, K); % MSC GPS aiding
    %         Dx = reshape(nav.Dx, N, K)';
    %         for k = 1:K
    %             nav.Xhat(k) = gpsins_correct_xhat(p, nav.Xhat(k), Dx(k,:)');
    %         end
    %         nav.Dx(:) = 0;
    %         nav.xhat = nav.Xhat(K);
    %         %nav.dx.dx_minus = Dx(K,:)';
    %         gps_plot_Xhat(p, nav.Xhat, nav.dx_buff, K, 200);
    %
    % %         %reset buffer
    % %         count = 0;
    % %         nav.Pxx=[];
    %     end
    
    
    pos_ecef = p.ecef_p_b + p.R_ned2ecef*nav.xhat.r_tb_t;
    [H, R, dy] = gpsins_gps_output(nav.param, pos_ecef, data.gps(ii_gps), 1); % always use DGPS not RTK
    [nav.dx, out] = gpsins_correct_dx(nav.dx, H, R, dy, out); % GPS aiding
    
    % correct navigation state
    if norm(nav.dx.dx_minus) > 0
        [nav.xhat] = gpsins_correct_xhat(nav.param, nav.xhat, nav.dx.dx_minus);
        count=count+1; % increment the buffer count
        if nav.dx.applied && count<=K
            %nav.rpos_dgps(count) = LS_double_diff(data.gps(ii_gps), pos_ecef, 1);
            nav.Xhat(count) = nav.xhat;
            nav.dx_buff(count) = nav.dx;
            nav.Pxx = blkdiag(nav.Pxx,nav.dx.Pxx_minus);
        elseif ~nav.dx.applied % if not correction happened
            count = max(0, count-1); % revert count
        end
        nav.dx.dx_minus(:) = 0;
    end
    if count>=K % if buffer full
        if count>K
            nav.Xhat = [nav.Xhat(2:end);nav.xhat]; % push xhat buffer
            nav.dx_buff = [nav.dx_buff(2:end), nav.dx]; % push dx buffer
        end
        
        %gps_plot_Xhat(p, nav.Xhat, nav.dx_buff, K, 100);
        
        Pos_ecef = zeros(K*3,1);
        for n=1:K
            Pos_ecef(3*n-2:3*n) = p.ecef_p_b + p.R_ned2ecef*nav.Xhat(n).r_tb_t;
        end
        [H, R, dy] = msc_gps_ins_output(nav.param, Pos_ecef, data.gps(ii_gps-K+1:ii_gps), K);
        [nav] = msc_gps_correct_Dx(nav, H, R, dy, K); % MSC GPS aiding
        Dx = reshape(nav.Dx, N, K)';
        for k = 1:K
            nav.Xhat(k) = gpsins_correct_xhat(p, nav.Xhat(k), Dx(k,:)');
        end
        nav.Dx(:) = 0;
        nav.xhat = nav.Xhat(K);
        %nav.dx.dx_minus = Dx(K,:)';
        %gps_plot_Xhat(p, nav.Xhat, nav.dx_buff, K, 200);
        
        %         %reset buffer
        %         count = 0;
        %         nav.Pxx=[];
    end   
end