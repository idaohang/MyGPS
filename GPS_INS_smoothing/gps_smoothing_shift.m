% This function does smoothing as
% 0. Decide the dimension of integer vector Nhat
% 1. First include the integers Nhat as a real vector, run Gauss-Newton
%    the initial value is the least square (average) result from Xhat_k
% 2. With the initial value from 1, solve the integer
function [out, Xhat_k, Nhat_k] = gps_smoothing_shift(nav,data,Xhat_k,Phase_meas)
% figure out how many integers should be included in Nhat and
% give initial Nhat based on Xhat_k
%Phase_meas = validate_phase_meas(nav, data, Xhat_k);

global Dual_Freq;
dual_freq = Dual_Freq;

% Integer candidate number
int_num = 15;

% dimension of states
Ns = nav.param.X_STATES;

% number of states
K = nav.buff.kdim;

% dimension of correction
l = 3 + Phase_meas.num_ph_sv-1;

Xhat_c = Xhat_k;
Xhat_t = Xhat_c; % initialize temp Xhat_t

% initial Shift
Shift = zeros(3,1);
% init Nhat
Phase_meas_t = Phase_meas;

disp('--------SHIFT_smoothing Start-------')

itc = 1;
maxit =10;
tol = 0.1;
alpha = 1.d-4;
[rc, A_mils, gc, B_mils, Ndim] = gps_smth_res_J_shift(nav, data, Xhat_t, Phase_meas_t);
% rc(1:3) = 0;
% A_mils(1:3,:) = eye(3);

% add virtual measurements from positions of Xhat_c
% Pos_Cov = 1; % 0.25
% SigmaV = chol(inv(Pos_Cov));
% HV = SigmaV*eye(3);
% rV = SigmaV*zeros(3,1);

% rc = [rc; rV]; % last term K*rV
% A_mils = [A_mils;HV];
% B_mils = [B_mils;zeros(3,Ndim*(1+dual_freq))];
n_rc0 = norm(rc)

%
%plot_Xhat(nav, Xhat_t, (K+1), itc);
%pause();

% A_mils = Jc(:,1:3);
% B_mils = Jc(:,4:end);
[Shift, N_hat] = mils( A_mils, B_mils, rc, int_num);

%% Solution selection
% kicked = [];
% list = [1:int_num];
% for i = 1:int_num
%     if norm(Shift(:,i))>2.6
%         kicked = [kicked;i];
%     end
% end
% for i = 1: length(kicked)   
%     list = list( list~=kicked(i) );
% end
% int_num_temp = int_num - length(kicked);
% Shift = Shift(:,list);
% N_hat = N_hat(:,list);
% if norm( Shift(:,1) - Shift(:,2) ) > norm( Shift(:,1) - Shift(:,3) )
%     N_hat(:,2) = N_hat(:,3);
%     Shift(:,2) = Shift(:,3);
%     display('Kick middle solution!!')
% end
% if int_num_temp>=2
%     if norm( N_hat(:,1) ) > norm( N_hat(:,2) ) && norm(Shift(:,1))>norm(Shift(:,2))
%         N_hat(:,1) = N_hat(:,2);
%         Shift(:,1) = Shift(:,2);
%         display('Pick alternative solution!!')
%     end
% end
%%
flag = 1;
%while( flag )
    % correct Xhat
    for k = 1:size(Xhat_c,1)
        Xhat_t(k) = gpsins_correct_xhat_shift(nav.param, Xhat_c(k), Shift(:,1));
    end
    % correct Nhat
    Phase_meas_t.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + N_hat(1:Ndim,1);
    if dual_freq
        Phase_meas_t.Nhat_wd(2:end) = Phase_meas_t.Nhat_wd(2:end) + N_hat(Ndim+1:end,1);
    end
    
    itc = itc + 1;
    %plot_Xhat(nav, Xhat_t, (K+1), itc);
    %pause();
    [rc, A_mils, gc, B_mils] = gps_smth_res_J_shift(nav, data, Xhat_t, Phase_meas_t);
    %     rc(1:3) = 0;
    %     A_mils(1:3,:) = eye(3);
%     rV = SigmaV*Shift(:,1);
%     rc = [rc; rV]; % last term K*rV
%     A_mils = [A_mils;HV];
%     B_mils = [B_mils;zeros(3,Ndim*(1+dual_freq))];
    norm(rc)
    %     A_mils = Jc(:,1:3);
    %     B_mils = Jc(:,4:end);
    [Shift, N_hat] = mils( A_mils, B_mils, rc, int_num);
    kicked = [];
    list = [1:int_num];
    for i = 1:int_num
        if norm(Shift(:,i))>2
            kicked = [kicked;i];
        end
    end
    for i = 1: length(kicked)
        list = list( list~=kicked(i) );
    end
    int_num_temp = int_num - length(kicked);
    Shift = Shift(:,list);
    N_hat = N_hat(:,list);
    
%     if norm( Shift(:,1) )<= 0.05 || norm( Shift(:,2) )<= 0.05 || norm( Shift(:,3) )<= 0.05
%         flag = 0;
%     else
%         display('One MILS is not enough!!')
%         pause();
%         for k = 1:size(Xhat_c,1)
%             Xhat_t(k) = gpsins_correct_xhat_shift(nav.param, Xhat_c(k), Shift(:,1));
%         end
%         % correct Nhat
%         Phase_meas_t.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + N_hat(:,1);
%     end
    
    %     if norm(N_hat(:,1))==0
    %         % Fix integer, optimizae X_hat
    %         [rc, Jc, gc, gps_res] = gps_smth_res_J_int(nav,data, Xhat_t, Phase_meas_t);
    %     end
%end

Xhat_k = Xhat_t;
Nhat_k = Phase_meas_t.Nhat; %
%plot_Xhat(nav, Xhat_k, (K+1), itc);

disp('--------SHIFT_smoothing End-------')

%% Compute integer
[rpos_ecef, delta_t] = LS_single_diff(data.gps(1), nav.param.ecef_p_b, 0);
rpos_tp=nav.param.R_ned2ecef' *( rpos_ecef - nav.param.ecef_p_b);
if ( norm(rpos_tp - Xhat_k(2,1).r_tb_t)>0.1)
    display('positioning error is too large!!')
end
data_struct = data.gps(1,1);
prn = data_struct.prnlist(1);
N_common = round( (  data_struct.Sat_state( 1, prn ).sd_phase_l1 * nav.param.wave_l1 - norm( rpos_ecef - data_struct.Sat_state( 1, prn ).sv_pos_ecef' ) ...
    - delta_t )/nav.param.wave_l1 );
delta_N = N_common - Nhat_k(1);
Nhat_k = Nhat_k + delta_N;

%%


for ii = 1:(K+1)
    % log data
    out.t(ii,:) = Xhat_k(ii).t;
    out.r_tb_t(ii,:) = (Xhat_k(ii).r_tb_t)';
    out.v_tb_t(ii,:) = (Xhat_k(ii).v_tb_t)';
    out.v_tb_b(ii,:) = (Xhat_k(ii).R_b2t'*Xhat_k(ii).v_tb_t)';
    out.R_b2t(:,:,ii)  = Xhat_k(ii).R_b2t;
    out.ba_b(ii,:) = Xhat_k(ii).ba_b';
    out.bg_b(ii,:) = Xhat_k(ii).bg_b';
    %out.Pxx(ii,:) = diag(nav.dx.Pxx_minus);
    xt = Xhat_k(ii).R_b2t*[1;0;0];
    out.alphahat(ii,:) = atan2(xt(2),xt(1));
    %out.res = gps_res;
end
%
% h = figure(1);
% subplot(3,1,1);
% plot(out.t, out.r_tb_t(:,1)-nav.param.base_pos(1),'--*r');
% title(['NED position'] );
% hold on;
% grid on;
% ylabel('x (m)');
% subplot(3,1,2);
% plot(out.t, out.r_tb_t(:,2)-nav.param.base_pos(2),'--*r');
% hold on;
% grid on;
% ylabel('y (m)');
% subplot(3,1,3);
% plot(out.t, out.r_tb_t(:,3)-nav.param.base_pos(3),'--*r');
% hold on;
% grid on;
% ylabel('z (m)');

end