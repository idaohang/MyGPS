function [H, R, residual] = msc_gps_output(p, rpos_ecef, meas_buff, K)
 k=K;
 meas_temp = meas_buff(k);
% rpos_dgps = LS_double_diff(meas_temp, rpos_ecef(3*k-2:3*k), 1);
% rpos_tp = p.R_ned2ecef' *( rpos_dgps - p.ecef_p_b) - p.true_pos
% %dis = rpos_dgps - rpos_ecef(3*k-2:3*k)
 [num_sat, Meas, R_meas, prn_dd] = gpsins_calc_meas_msc(meas_temp); % get the latest l1 phase in meters
% Meas = Meas*p.wave_l1;
H = zeros(num_sat*K,K*3); % 3 is the state(position) size

residual = zeros(num_sat*K,1);
R = zeros(num_sat*K,1);
 Nhat = zeros(num_sat,1);
%V = ones(num_sat*K,1);
V = zeros(num_sat*K,num_sat);
for i=1:num_sat
    V(K*(i-1)+1:K*i,i) = p.wave_l1*ones(K,1);
end
%V =ones(K,1);
[Q, R_temp] = qr(V);
%A = Q(:,num_sat-2:num_sat*K);
A = Q(:,num_sat+1:num_sat*K);
%A = Q(:,2:K);
%A=1;
for i=1:num_sat
    prn = prn_dd(i,:);
%      Nhat(i) = round((Meas(i)-norm(rpos_ecef(3*k-2:3*k)'-meas_temp.Sat_state(prn(1)).sv_pos_ecef)...
%     +norm(rpos_ecef(3*k-2:3*k)'-meas_temp.Sat_state(prn(2)).sv_pos_ecef))/p.wave_l1);

    Nhat(i) = round((Meas(i)-norm(rpos_dgps'-meas_temp.Sat_state(prn(1)).sv_pos_ecef)...
    +norm(rpos_dgps'-meas_temp.Sat_state(prn(2)).sv_pos_ecef))/p.wave_l1);
    % for debug only
     Nhat(i) = -(meas_temp.Sat_state(prn(1)).N_l1-meas_temp.Sat_state(prn(2)).N_l1);

% for debug only
%     Nhat(i)
%     round((Meas(i)-meas_temp.Sat_state(prn(1)).ph_l1_rng...
%     +meas_temp.Sat_state(prn(2)).ph_l1_rng)/p.wave_l1)

%     Nhat(i) = round((Meas(i)-meas_temp.Sat_state(prn(1)).sd_code_l1....
%     +meas_temp.Sat_state(prn(2)).sd_code_l1)/p.wave_l1);
end


for j=1:K
    meas_temp = meas_buff(j); % for each measurement instant, from latest to oldest    
    [num_sat, sd_phase_l1, R_meas, prn_dd] = gpsins_calc_meas_msc(meas_temp); % get the l1 phase in meters
    Meas = (sd_phase_l1 )*p.wave_l1; % compensate integer ambiguity
    for i=1:num_sat
        [Hij, PRangei] = gpsins_sat_output(rpos_ecef(3*j-2:3*j), prn_dd(i,2), prn_dd(i,1), meas_temp);
        H((i-1)*K+j,3*j-2:3*j) = Hij;%
        residual((i-1)*K+j) = (Meas(i) - PRangei);
        R((i-1)*K+j) = R_meas(i,i);
        
%         H((j-1)*num_sat+i,(3*j-2):3*j)= Hij*p.R_ned2ecef;
%         residual((j-1)*num_sat+i) = Meas(i) - PRangei;
%         R((j-1)*num_sat+i,(j-1)*num_sat+i) = R_meas(i,i);        
    end
end
H = A'*H;
residual = A'*residual;
R=A'*diag(R)*A;

% H0 = [];
% r0 = [];
% R0 = [];
% for i = 1:num_sat
%     H_temp = A'* H(K*(i-1)+1:K*i,:);
%     r_temp = A'* residual(K*(i-1)+1:K*i);
%     R_temp = A'*diag(R((K*(i-1)+1:K*i)))*A;
%     H0 = [H0;H_temp];
%     r0 = [r0;r_temp];
%     R0 = blkdiag(R0,R_temp);
% end
% H = H0;
% residual = r0;
% R = R0;
end



% output residuals for analysis
% 	if ~isempty(out)
% 		if isfield(out, 'residual') & idx <= length(out.residual) & isfield(out.residual{idx}, 't')
% 			i = length(out.residual{idx}.t) + 1;
% 		else
% 			i = 1;
% 		end
% 		out.residual{idx}.t(i,:) = dx.t';
% 		out.residual{idx}.HPHtR(i,:) = diag(HPHtR)';
% 		out.residual{idx}.dy(i,:) = dy';
% 		out.residual{idx}.dx_plus(i,:) = dx_plus';
% 	end