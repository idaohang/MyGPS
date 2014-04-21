% This function gives absolute phase measurement of i-th satellite on the
% prnlist, in L1 frequency
% p is the parameter
% rpos_ecef are the positions 
% meas_buff is the buffer of gps data
% K is the window length (seconds)
% i is the index of satellite on the prn list
% Phase_meas the prn list and sat number

% Output is Jacobian of position, R matrix for noise, dy, R1 is Jacobian
% for Ni
function [H, R, residual, R1] = abs_phase_output(p, rpos_ecef, meas_buff, i, K, Phase_meas)
global Dual_Freq;
dual_freq = Dual_Freq;
% init output
H = zeros(K, 3); % 3 is for position dimension
residual = zeros(K,1);
R = zeros(K,1);

if dual_freq
    residual_wd = zeros(K,1);
end

V = ones(K,1);
[Q1, R_1] = qr(V);
R1 = R_1(1,:)*p.wave_l1;
R_wd = R_1(1,:)*p.wave_wd;
B = Q1(:,1);

for j=1:K
    meas_temp = meas_buff(j); % for each measurement instant, from latest to oldest    
    [num_sat, dd_phase, R_meas, prn_dd] = msc_gps_calc_meas(meas_temp, Phase_meas.prn_phase); % get the l1 phase in cycles
    if num_sat ~= (Phase_meas.num_ph_sv - 1)
        error('Number of Valid Satellite does not match!!')
    end
    Meas = dd_phase(1:num_sat)*p.wave_l1; 
    
    if Phase_meas.prn_phase(i) ~=  prn_dd(i-1,1) || Phase_meas.prn_phase(1) ~= prn_dd(i-1,2)
        error('Common satellite or picked satellite does not match!!')
    end
    [Hij, PRangei] = gpsins_sat_output(rpos_ecef(3*j-2:3*j), prn_dd(i-1,2), prn_dd(i-1,1), meas_temp);
    H(j,:) = Hij*p.R_ned2ecef;
    residual(j) = Meas(i-1) - PRangei;
    R(j) = R_meas(i-1,i-1);
    
    if dual_freq
        Meas_wd = ( dd_phase(1:num_sat) - dd_phase(num_sat+1:end) )*p.wave_wd;
        residual_wd(j) = Meas_wd(i-1) - PRangei;
    end
end

residual = B'*residual -  R1*( Phase_meas.Nhat(i) - Phase_meas.Nhat(1) );
if dual_freq
    residual_wd = B'*residual_wd -  R_wd*( Phase_meas.Nhat_wd(i) - Phase_meas.Nhat_wd(1) );
    H = [ B'*H; B'*H];
    residual = [ residual; residual_wd ];
    R = [ B'*diag(R)*B; B'*diag(5.7^2*R)*B ];  % 5.7 is from JAF2008 below eqn.8.106
    R1 = [R1;R_wd];
else
    H = [B'*H];    
    R=B'*diag(R)*B;
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
