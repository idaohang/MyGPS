clc;
clear;
%% Options
Start_time = 17;
End_time = 97;

%% Load data
load GPS_buff_log22.dat;
buff = read_GPS_buff_log(GPS_buff_log22, Start_time, End_time);
%% Simulation setup
% simulation parameters
p = gpsins_param();

K = length(buff);  % number of states

rpos_ecef = zeros(3*K,1); % init rover state
true_pos = zeros(3*K,1);
error = zeros(3*K,1);
for i = 1:K
    ic = randomcov(100, eye(3), [0,0,0]); % generate 1000 random vector with covariance and mean [0,0,0] 
    error(3*i-2:3*i) = 10*ic(1,:);
    rpos_ecef(3*i-2:3*i) = p.ecef_p_b+p.R_ned2ecef*(p.true_pos+error(3*i-2:3*i));
    true_pos(3*i-2:3*i) = p.ecef_p_b+p.R_ned2ecef*p.true_pos;
end

%% Least Square
recur_count = 0;
while(recur_count==0 || norm(delta)> 0.000001)
    for i = 1:K
        dis(3*i-2:3*i) = p.R_ned2ecef'*(rpos_ecef(3*i-2:3*i)-true_pos(3*i-2:3*i));
    end
    fig = 1;
    h = figure(fig);
    subplot(3,1,1);
    plot(1:K, dis(1:3:3*K),'*');
    title(['NED position'] );
    grid on;
    ylabel('x (m)');
    subplot(3,1,2);
    plot(1:K, dis(2:3:3*K),'*');
    grid on;
    ylabel('y (m)');
    subplot(3,1,3);
    plot(1:K, dis(3:3:3*K),'*');
    grid on;
    ylabel('z (m)');
    
    [H, R, dy] = msc_gps_output(p, rpos_ecef, buff, K);
    norm(dy)
    Cov_meas = diag(1./diag(R));
    delta = (H'*Cov_meas*H)\(H'*Cov_meas*dy);
%     for i = 1:K
%         delta(3*i-2:3*i) = p.R_ned2ecef*delta(3*i-2:3*i); % transform to in ECEF frame
%     end
    rpos_ecef = rpos_ecef + delta;
    recur_count=recur_count+1;
end