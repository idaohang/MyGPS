% This function does smoothing as
% 0. Decide the dimension of integer vector Nhat
% 1. First include the integers Nhat as a real vector, run Gauss-Newton
%    the initial value is the least square (average) result from Xhat_k
% 2. With the initial value from 1, solve the integer
function [out, Xhat_k, Nhat_k] = gps_smoothing_int(nav,data,Xhat_k,Phase_meas)
% figure out how many integers should be included in Nhat and
% give initial Nhat based on Xhat_k
%Phase_meas = validate_phase_meas(nav, data,Xhat_k);

global Dual_Freq;
dual_freq = Dual_Freq;

% dimension of states
Ns = nav.param.X_STATES;

% number of states
K = nav.buff.kdim;

% dimension of correction
l = Ns*(K+1) + Phase_meas.num_ph_sv-1;

% initial Xhat
Xhat_c = Xhat_k;
Xhat_t = Xhat_c; % initialize temp Xhat_t
% init Nhat
Phase_meas_t = Phase_meas;
num_sat= (Phase_meas.num_ph_sv);

itc = 1;
maxit =10;
tol = 0.01;
alpha = 1.d-4;
[rc, Jc, gc, gps_res] = gps_smth_res_J_int(nav,data, Xhat_t, Phase_meas_t);
n_rc0 = norm(rc)
delta = 1;
%
% plot_Xhat(nav, Xhat_t, (K+1), itc);
% pause();

while(norm(delta)>tol && itc<=maxit)
    lambda = 1;
    % Gauss-Newton
    [Q, R] = qr(Jc);
    de = -Q'*rc;
    d = de(1:l, :);
    R = R(1:l, :);
    dc = R\d;
    dgn = lambda*dc;
    %n_dgn = norm(dgn);
    delta = reshape(dgn(1:Ns*(K+1)), Ns, (K+1))';
    % correct Xhat
    for k = 1:size(Xhat_c,1)
        Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
    end
    % correct Nhat
    Phase_meas_t.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + dgn(Ns*(K+1)+1:Ns*(K+1)+num_sat-1);
    if dual_freq
        Phase_meas_t.Nhat_wd(2:end) = Phase_meas_t.Nhat_wd(2:end) + dgn(Ns*(K+1)+num_sat:end);
    end
    
    %     n_rc0 = norm(rc)
    %     itc
    %     plot_Xhat(nav, Xhat_t, (K+1), itc+1);
    %     pause()
    
    [rt, Jt, gt] = gps_smth_res_J_int(nav, data, Xhat_t, Phase_meas_t);
    iarm=0;
    itc=itc+1;
    %
    % Goal for sufficient decrease
    %
    rgoal = norm(rc) - alpha*lambda*(gc'*dc); % dgn=lambda*dc
    Armijo = 0;
    while(norm(rt)>rgoal)
        iarm=iarm+1;
        lambda=lambda/2;
        rgoal = norm(rc) - alpha*lambda*(gc'*dc);
        dgn = lambda*dc;
        delta = reshape(dgn(1:Ns*(K+1)), Ns, (K+1))';
        for k = 1:size(Xhat_c,1)
            Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
        end
        % correct Nhat
        Phase_meas_t.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + dgn(Ns*(K+1)+1:Ns*(K+1)+num_sat-1);
        if dual_freq
            Phase_meas_t.Nhat_wd(2:end) = Phase_meas_t.Nhat_wd(2:end) + dgn(Ns*(K+1)+num_sat:end);
        end
        
        [rt, Jt, gt] = gps_smth_res_J_int(nav, data, Xhat_t, Phase_meas_t);
        if(iarm > 10)
            disp(' Armijo error in Gauss-Newton')
            Armijo = 1;
            break;
        end
    end
    if Armijo
        break;
    else
        Xhat_c = Xhat_t;
        Phase_meas = Phase_meas_t;
        [rc, Jc, gc, gps_res] = gps_smth_res_J_int(nav, data, Xhat_c, Phase_meas);
        norm(rc)
    end
end
Xhat_k = Xhat_c;
Nhat_k = Phase_meas.Nhat;
plot_Xhat(nav, Xhat_k, (K+1), itc);

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
    out.res = gps_res;
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