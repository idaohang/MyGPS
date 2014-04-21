%%
% 0119_1042 Veh_data_0613_1055
% 20140124_1728 Ahn's loop
% 20140124_1848 600s open sky
% 20140311_1055 1000s open skyh

clc;
clear;
%close all;

global DGPS;
global Dual_Freq;
global Stationary;
global Elev_Mask;
global Sat_Num;
DGPS = 1;      % only use single differenced code meas
Dual_Freq = 1; 
Stationary = 0;
Elev_Mask = 5;
Sat_Num = 15;
global block_list;
block_list = [];%[32, 16];

% for stationary
true = [0.025, -5.931, 0.051];
%% Add routine path
addpath('../');
addpath('../MILES');
addpath('../Data_set')
addpath('../util')

global imu_data;
global gps_data;
if ~Stationary
%     imu_data = load('../Data_set/IMU_buff_log_veh.dat');
%     gps_data = load('../Data_set/GPS_buff_log_veh.dat');
    %IMU_meas_log20140124_1848;
    IMU_meas_log20140119_1042;
    imu_data = imu_meas;
    %gps_data = load('GPS_buff_log20140124_1848.dat');
    gps_data = load('GPS_buff_log20140119_1042.dat');
    gps_data(:,1) = gps_data(:,2);
    load('../Data_set/python_smooth.mat'); % smoothing result in state_out
    %load('../Data_set/mat_smoothing_600s.mat');
    %load('Ahns_loop.mat'); % smoothing result in state_out
    %state_out = state_out;
else
    imu_data = load('../Data_set/IMU_buff_log_UCRStationary_1000s.dat');
    gps_data = load('../Data_set/GPS_buff_log_UCRStationary_1000s.dat');
    %gps_data = load('GPS_buff_log20140207_0728.dat');
end

[K,~] = size(gps_data);
skipp = 5;     % skip the first xx seconds
buff_size = 5;
% overlap = 0.2;
% gap = buff_size*overlap;
% test_num = round((K - skipp)/(buff_size*overlap));
gap = 1;
check_list = skipp:gap:K;
%check_list = sort([42]);
test_num = length(check_list);
Result = zeros(test_num, 30); % 10
p = gpsins_param();

for jj = 1:test_num  %skipp:buff_size*overlap:K   
    k = check_list(jj);
    start_time = round( gps_data(k,1) );
    if (start_time+buff_size)>gps_data(end,1);
        display('Reached the end.');
        break;
    end
    
    if ~Stationary
        % trim data for plotting
        i = 1;
        while ( state_out(i,1) < start_time )
            i=i+1;
        end
        j=i+1;
        while ( state_out(j,1)<start_time+buff_size-0.5)
            j=j+1;
        end
        state_out_k = state_out(i:j,:);
%         [l,~] = size(state_out_k);
%         for i=1:l
%             state_out_k(i,2:4) = p.R_ned2ecef'*(state_out_k(i,2:4)' - p.ecef_p_b);
%         end
    else
        state_out_k = [];
    end
    
    [N1, N2, N3, n, X1, X2, X3, out1, out2, out3, out0, mils_res] = gps_smth_ins( start_time , buff_size, state_out_k );
    
    m = size(state_out_k);
    if m
        [error0, var0, std0] = check_result_error(out0, state_out_k); % int-free
        [error1, var1, std1] = check_result_error(out1, state_out_k, 1); % shift
        [eee, ~] = check_result_error(out1, out2, 1) % shift
        [error2, var2, std2] = check_result_error(out2, state_out_k, 1); % mils
        [error3, var3, std3] = check_result_error(out3, state_out_k); % floating
        [mm,nn] = size(mils_res);
        Mils_Result((jj-1)*buff_size+1:jj*buff_size,1:nn) = mils_res;
    else       
        out11 = out1;
        one_vec = ones(length(out11),1);
        ground_truth = [true(1)*one_vec, true(2)*one_vec, true(3)*one_vec];
        error1 = mean( ground_truth - out11(:,2:4));
        var1 = var( ground_truth - out11(:,2:4));
        std1 = std( ground_truth - out11(:,2:4));
        
        out12 = out2;
        error2 = mean( ground_truth - out12(:,2:4));
        var2 = var( ground_truth - out12(:,2:4));
        std2 = std( ground_truth - out12(:,2:4));
        
        out13 = out3;
        error3 = mean( ground_truth - out13(:,2:4));
        var3 = var( ground_truth - out13(:,2:4));
        std3 = std( ground_truth - out13(:,2:4));
        
        out10 = out0;
        error0 = mean( ground_truth - out10(:,2:4));
        var0 = var( ground_truth - out10(:,2:4));
        std0 = std( ground_truth - out10(:,2:4));
        
        [mm,nn] = size(mils_res);
        Mils_Result((k-skipp)*buff_size+1:(k-skipp+1)*buff_size,1:nn) = mils_res;
    end
    display(['Pos int-free mean is ', num2str(error0)] );
    display(['Pos shift mean is ', num2str(error1)] );
    display(['Pos mils mean is ', num2str(error2)] );
    display(['Pos floating mean is ', num2str(error3)] );
    display(['Pos int-free variance is ', num2str(var0)] );
    display(['Pos shift variance is ', num2str(var1)] );
    display(['Pos mils variance is ', num2str(var2)] );
    display(['Pos floating variance is ', num2str(var3)] );
    display(['Pos int-free std_dev is ', num2str(std0)] );
    display(['Pos shift std_dev is ', num2str(std1)] );
    display(['Pos mils std_dev is ', num2str(std2)] );
    display(['Pos floating std_dev is ', num2str(std3)] );
    
    if length(N1) > 6
        Result(jj,:) = [start_time, start_time+buff_size, ...
            n, error0, std0, sqrt(norm(N1)^2/length(N1)), error1, std1, sqrt(norm(N2)^2/length(N2)), ...
            error2, std2, sqrt(norm(N3)^2/length(N2)), error3, std3 ];
    end
    %Result((k-skipp)/(buff_size*overlap)+1,:) = [start_time, start_time+buff_size, norm(N1), n, error1, std1];

    
    % for stationary
%     true = [0.028, -5.92, 0.05];
%     Pos_tp1 = zeros( (buff_size+1), 3 );
%     Pos_tp2 = zeros( (buff_size+1), 3 );
%     for k=1:buff_size+1
%         Pos_tp1(k,:) = X1(k).r_tb_t';
%         Pos_tp2(k,:) = X2(k).r_tb_t';
%     end
%     error1 = mean(out1.r_tb_t - true);
%     error2 = mean(out2.r_tb_t - true);
%     Result(i-skipp+1,:) = [ norm(N1), norm(N2), n, error1, error2, start_time, start_time+buff_size ];
end