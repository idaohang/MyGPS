% This function does smoothing as
% 0. Decide the dimension of integer vector Nhat
% 1. First include the integers Nhat as a real vector, run Gauss-Newton
%    the initial value is the least square (average) result from Xhat_k
% 2. With the initial value from 1, solve the integer
function [out, Xhat_k, Nhat_k, gps_res] = gps_smoothing_miles(nav,data,Xhat_k,Phase_meas)
% figure out how many integers should be included in Nhat and
% give initial Nhat based on Xhat_k
%Phase_meas = validate_phase_meas(nav, data, Xhat_k);

global Dual_Freq;
dual_freq = Dual_Freq;

% dimension of states
Ns = nav.param.X_STATES;

% number of states
K = nav.buff.kdim;

% dimension of correction
l = Ns*(K+1) + Phase_meas.num_ph_sv-1;

% initial Xhat
Xhat_o = Xhat_k; % save the original estimate
Xhat_c = Xhat_k;
Xhat_t = Xhat_c; % initialize temp Xhat_t
% init Nhat
Phase_meas_0 = Phase_meas;
Phase_meas_t = Phase_meas;
num_sat= (Phase_meas.num_ph_sv);

itc = 1;
[rc, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
n_rc0 = norm(rc)
delta = 1;

% 
%plot_Xhat(nav, Xhat_t, (K+1), itc);
%pause();

disp('--------MILS_smoothing Start-------')

A_mils = Jc(:,1:Ns*(K+1));
B_mils = Jc(:,Ns*(K+1)+1:end);
[X_Hat, N_hat] = mils( A_mils, B_mils, rc, 10);

delta = reshape(X_Hat(:,1), Ns, (K+1))';

while( norm(N_hat(:,1))>0.8 && itc<10)   
    % correct Xhat
    for k = 1:size(Xhat_t,1)
        Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_t(k), delta(k,:)');
    end
    % correct Nhat
    Phase_meas_t.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + N_hat(1:num_sat-1,1);
    if dual_freq
        Phase_meas_t.Nhat_wd(2:end) = Phase_meas_t.Nhat_wd(2:end) + N_hat(end-num_sat+2:end,1);
    end
    
    itc = itc + 1
%     plot_Xhat(nav, Xhat_t, (K+1), itc);
%     pause();
    [rc, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
    norm(rc)
    A_mils = Jc(:,1:Ns*(K+1));
    B_mils = Jc(:,Ns*(K+1)+1:end);
    [X_Hat, N_hat] = mils( A_mils, B_mils, rc, 3);
    
    delta = reshape(X_Hat(:,1), Ns, (K+1))';
end
%plot_Xhat(nav, Xhat_t, (K+1), 3);


Xhat_c = Xhat_t;
% Then fix N, minimize X
itc = 1;
maxit = 30;
tol = 0.02;
alpha = 1.d-4;
[rc, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
norm(rc)
Jc = Jc(:,1:Ns*(K+1));
while (norm(delta)>tol && itc<=maxit)
    lambda = 1;
    % Gauss-Newton
    [Q, R] = qr(Jc);
    de = -Q'*rc;
    d = de(1:l, :);
    R = R(1:l, :);
    dc = R\d;
    %dc = -(Jc'*Jc)\(Jc'*rc);
    dgn = lambda*dc;
    %n_dgn = norm(dgn);
    delta = reshape(dgn, Ns, K+1)';
    
    for k = 1:size(Xhat_c,1)
        Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
    end
%     n_rc0 = norm(rc)
%     itc
    
    
%     plot_Xhat(nav, Xhat_t, K, itc+1);
%     pause()
    
    [rt, Jt, gt, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
    Jt = Jt(:,1:Ns*(K+1));
    iarm=0;
    itc=itc+1;
    %
    % Goal for sufficient decrease
    %    
    %rgoal = norm(rc) - alpha*lambda*(gc'*dc); % dgn=lambda*dc
    Armijo = 0;
    %while(norm(rt)>rgoal)
    while(norm(rt)>norm(rc))
        iarm=iarm+1;
        lambda=lambda/2;  % 1.2
        %rgoal = norm(rc) - alpha*lambda*(gc'*dc);
        dgn = lambda*dc;
        delta = reshape(dgn, Ns, K+1)';
        for k = 1:size(Xhat_c,1)
            Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
        end
        [rt, Jt, gt, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
        Jt = Jt(:,1:Ns*(K+1));
        if(iarm > 10)
            disp(' Armijo error in Gauss-Newton for MILS')
            Armijo = 1;
            break;
        end
    end
    if Armijo
        break;
    else
        Xhat_c = Xhat_t;
        [rc, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
        Jc = Jc(:,1:Ns*(K+1));
%         norm(rc)
    end        
end
disp(['dX = ', num2str(norm(delta))]);
disp(['iteration ', num2str(itc)]);
disp(['cost ', num2str(norm(rc))]);
disp('--------MILS_smoothing End-------')
% [rc, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);
% norm(rc)
% A_mils = Jc(:,1:Ns*(K+1));
% B_mils = Jc(:,Ns*(K+1)+1:end);
% [X_Hat, N_hat] = mils( A_mils, B_mils, rc, 3);

Xhat_k = Xhat_c;
Nhat_k = Phase_meas_t.Nhat;  %(2:end)
%plot_Xhat(nav, Xhat_k, (K+1), 3);
[rc_1, Jc, gc, gps_res] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t);

%% Compute second candidate for Ratio Tests
% A_mils = Jc(:,1:Ns*(K+1));
% B_mils = Jc(:,Ns*(K+1)+1:end);
% [X_Hat, N_hat] = mils( A_mils, B_mils, rc_1, 10);
% 
% delta = reshape(X_Hat(:,1), Ns, (K+1))';
% 
% % Prepare Second candidate for Ratio Test
% delta_2 = reshape(X_Hat(:,2), Ns, (K+1))';
% % Ratio Test
% Phase_meas_t_2 = Phase_meas_t;
% Xhat_temp2 = Xhat_t;
% for k = 1:size(Xhat_c,1)
%     Xhat_temp2(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta_2(k,:)');
% end
% Phase_meas_t_2.Nhat(2:end) = Phase_meas_t.Nhat(2:end) + N_hat(1:num_sat-1,2);
% if dual_freq
%     Phase_meas_t_2.Nhat_wd(2:end) = Phase_meas_t.Nhat_wd(2:end) + N_hat(end-num_sat+2:end,2);
% end
% 
% itc = 1;
% maxit =30;
% tol = 0.02;
% delta = 1;
% alpha = 1.d-4;
% [rc, Jc, gc, gps_res_2] = gps_smth_res_J_mils(nav,data, Xhat_temp2, Phase_meas_t_2);
% norm(rc)
% Jc = Jc(:,1:Ns*(K+1));
% while (norm(delta)>tol && itc<=maxit)
%     lambda = 1;
%     % Gauss-Newton
%     [Q, R] = qr(Jc);
%     de = -Q'*rc;
%     d = de(1:l, :);
%     R = R(1:l, :);
%     dc = R\d;
%     %dc = -(Jc'*Jc)\(Jc'*rc);
%     dgn = lambda*dc;
%     %n_dgn = norm(dgn);
%     delta = reshape(dgn, Ns, K+1)';
%     
%     for k = 1:size(Xhat_c,1)
%         Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_temp2(k), delta(k,:)');
%     end
%     n_rc0 = norm(rc)
%     itc
%     
%     
% %     plot_Xhat(nav, Xhat_t, K, itc+1);
% %     pause()
%     
%     [rt, Jt, gt, gps_res_2] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t_2);
%     Jt = Jt(:,1:Ns*(K+1));
%     iarm=0;
%     itc=itc+1;
%     %
%     % Goal for sufficient decrease
%     %    
%     %rgoal = norm(rc) - alpha*lambda*(gc'*dc); % dgn=lambda*dc
%     Armijo = 0;
%     %while(norm(rt)>rgoal)
%     while(norm(rt)>norm(rc))
%         iarm=iarm+1;
%         lambda=lambda/2;  % 1.2
%         %rgoal = norm(rc) - alpha*lambda*(gc'*dc);
%         dgn = lambda*dc;
%         delta = reshape(dgn, Ns, K+1)';
%         for k = 1:size(Xhat_c,1)
%             Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_temp2(k), delta(k,:)');
%         end
%         [rt, Jt, gt, gps_res_2] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t_2);
%         Jt = Jt(:,1:Ns*(K+1));
%         if(iarm > 10)
%             disp(' Armijo error in Gauss-Newton for MILS')
%             Armijo = 1;
%             break;
%         end
%     end
%     if Armijo
%         break;
%     else
%         Xhat_temp2 = Xhat_t;
%         [rc, Jc, gc, gps_res_2] = gps_smth_res_J_mils(nav,data, Xhat_t, Phase_meas_t_2);
%         Jc = Jc(:,1:Ns*(K+1));
%         norm(rc)
%     end        
% end
% Ratio = norm(rc)/norm(rc_1)

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