function [H, R, residual] = msc_gps_ins_output(p, rpos_ecef, meas_buff, K, prn_list)
global Dual_Freq;
dual_freq = Dual_Freq;
k=1;
meas_temp = meas_buff(k);
% rpos_dgps = LS_double_diff(meas_temp, rpos_ecef(3*k-2:3*k), 1);
% rpos_tp = p.R_ned2ecef' *( rpos_dgps - p.ecef_p_b);
[num_sat, Meas, R_meas, prn_dd] = msc_gps_calc_meas(meas_temp, prn_list); % get sat num to preallocate
%Meas = Meas*p.wave_l1;
H = zeros(num_sat*K,K*p.X_STATES); 

if ~dual_freq
    res_temp = zeros(num_sat*K,1);
else
    res_temp = zeros(num_sat*K,2);
end
R = zeros(num_sat*K,1);

%Nhat = zeros(num_sat,1);
%V = ones(K,1);
V = zeros(num_sat*K,num_sat);
for i=1:num_sat
    V(K*(i-1)+1:K*i,i) = ones(K,1);
end
[Q, R_temp] = qr(V);
%A = Q(:,2:K);
A = Q(:,num_sat+1:num_sat*K);
for i=1:num_sat
     prn = prn_dd(i,:);
%     if unlocked(i)
%         Nhat(i) = round((Meas(i)-norm(rpos_dgps'-meas_temp.Sat_state(prn(1)).sv_pos_ecef)...
%             +norm(rpos_dgps'-meas_temp.Sat_state(prn(2)).sv_pos_ecef))/p.wave_l1);
%         unlocked(i) = 0;
%     else
%         Nhat(i) = round((Meas(i)-norm(rpos_ecef(3*k-2:3*k)'-meas_temp.Sat_state(prn(1)).sv_pos_ecef)...
%             +norm(rpos_ecef(3*k-2:3*k)'-meas_temp.Sat_state(prn(2)).sv_pos_ecef))/p.wave_l1);
%     end
    % for debug only
    %Nhat(i) = -(meas_temp.Sat_state(prn(1)).N_l1-meas_temp.Sat_state(prn(2)).N_l1);
    
% for debug only
%     Nhat(i)
%     round((Meas(i)-meas_temp.Sat_state(prn(1)).ph_l1_rng...
%     +meas_temp.Sat_state(prn(2)).ph_l1_rng)/p.wave_l1)

%     Nhat(i) = round((Meas(i)-meas_temp.Sat_state(prn(1)).sd_code_l1....
%     +meas_temp.Sat_state(prn(2)).sd_code_l1)/p.wave_l1);
end
%Ndiff = (Nhat - Nhat1')

for j=1:K
    meas_temp = meas_buff(j); % for each measurement instant, from latest to oldest    
    [num_sat, sd_phase_meas, R_meas, prn_dd] = msc_gps_calc_meas(meas_temp, prn_list); % get the l1 phase in meters
    Meas = sd_phase_meas(1:num_sat)*p.wave_l1; % phase meas in meters
    if dual_freq
        Meas = [ Meas; ( sd_phase_meas(1:num_sat) - sd_phase_meas(1+num_sat:end) )*p.wave_wd];
    end
    for i=1:num_sat
        [Hij, PRangei] = gpsins_sat_output(rpos_ecef(3*j-2:3*j), prn_dd(i,2), prn_dd(i,1), meas_temp);
        H((i-1)*K+j,p.X_STATES*(j-1)+p.X_POS) = Hij*p.R_ned2ecef;
        res_temp((i-1)*K+j,1) = Meas(i) - PRangei ;%- Nhat(i)*p.wave_l1;
        R((i-1)*K+j) = R_meas(i,i);
        if dual_freq
            res_temp((i-1)*K+j,2) = Meas(i+num_sat) - PRangei;
        end
        
%         H((j-1)*num_sat+i,(3*j-2):3*j)= Hij*p.R_ned2ecef;
%         residual((j-1)*num_sat+i) = Meas(i) - PRangei;
%         R((j-1)*num_sat+i,(j-1)*num_sat+i) = R_meas(i,i);        
    end
end

residual = A'*res_temp(:,1);
if dual_freq
    H = [ A'*H; A'*H];
    residual = [ residual; A'*res_temp(:,2) ];
    R = blkdiag( A'*diag(R)*A, A'*diag(5.7^2*R)*A);  % 5.7 is from JAF2008 below eqn.8.106
else
    H = A'*H;
    R = A'*diag(R)*A;
end

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
