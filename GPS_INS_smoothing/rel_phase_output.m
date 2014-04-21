% p is the parameter
% rpos_ecef are the positions 
% meas_buff is the buffer of gps data
% K is the window length (seconds)
% Phase_meas the prn list and sat number

% Output is Jacobian of position, R matrix for noise, dy, R1 is Jacobian
% for Ni
function [H, R, residual, R1] = rel_phase_output(p, rpos_ecef, meas_buff, i, K, Phase_meas)

% init output
H = zeros(K, 3); % 3 is for position dimension
residual = zeros(K,1);
residual1 = zeros(K,1);
R = zeros(K,1);

V = ones(K,1);
[Q1, R1] = qr(V);
R1 = R1(1,:)*p.wave_l1;
A = Q1(:,2:end);

for j=1:K
    meas_temp = meas_buff(j); % for each measurement instant, from latest to oldest    
    [num_sat, sd_phase_l1, R_meas, prn_dd] = msc_gps_calc_meas(meas_temp); % get the l1 phase in meters
    if num_sat ~= (Phase_meas.num_ph_sv - 1)
        error('Number of Valid Satellite does not match!!')
    end
    Meas = (sd_phase_l1)*p.wave_l1; 
    
    if Phase_meas.prn_phase(i) ~=  prn_dd(i-1,1) || Phase_meas.prn_phase(1) ~= prn_dd(i-1,2)
        error('Common satellite or picked satellite does not match!!')
    end
    [Hij, PRangei] = gpsins_sat_output(rpos_ecef(3*j-2:3*j), prn_dd(i-1,2), prn_dd(i-1,1), meas_temp);
    H(j,:) = Hij*p.R_ned2ecef;
    residual(j) = Meas(i-1) - PRangei;
    residual1(j) = Meas(i-1) - PRangei - p.wave_l1*( Phase_meas.Nhat(i) - Phase_meas.Nhat(1) );
    R(j) = R_meas(i-1,i-1);
end

H = [A'*H];
residual = A'*residual;
residual1 = A'*residual1;
R=A'*diag(R)*A;

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
