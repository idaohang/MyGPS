function [HH, R_meas, residual] = gpsins_gps_output_int(p, rpos_ecef, meas_temp, Phase_meas)
global Dual_Freq;
dual_freq = Dual_Freq;
[Meas, R_meas] = gpsins_calc_meas_int(meas_temp, Phase_meas);
num_meas = Phase_meas.num_ph_sv - 1;
if dual_freq
    Meas = [Meas(1:num_meas,1).*p.wave_l1; (Meas(1:num_meas,1) - Meas(num_meas+1:end,1)).*p.wave_wd];
else
    Meas = Meas(1:num_meas,1).*p.wave_l1;   % converstion from cycles to meters
end

svpos_ecef = zeros(3,1);
H0 = zeros(1,3);
H = zeros(num_meas,3);
Rng_comp = zeros(num_meas,1);
PRange = zeros(num_meas*(dual_freq+1),1);
HH = zeros(num_meas,p.X_STATES);

prn0 = Phase_meas.prn_phase(1);
svpos_ecef = meas_temp.Sat_state(prn0).sv_pos_ecef;
H0(1) = svpos_ecef(1) - rpos_ecef(1);
H0(2) = svpos_ecef(2) - rpos_ecef(2);
H0(3) = svpos_ecef(3) - rpos_ecef(3);
Rng_comp0 = norm( H0 );
H0(1) = H0(1)/(-Rng_comp0);
H0(2) = H0(2)/(-Rng_comp0);
H0(3) = H0(3)/(-Rng_comp0);

for i=1:num_meas    
    prni = Phase_meas.prn_phase(i+1);
    svpos_ecef = meas_temp.Sat_state(prni).sv_pos_ecef;
    H(i,1) = svpos_ecef(1) - rpos_ecef(1);
    H(i,2) = svpos_ecef(2) - rpos_ecef(2);
    H(i,3) = svpos_ecef(3) - rpos_ecef(3);
    Rng_comp(i) = norm( H(i,:) );
    H(i,1) = H(i,1)/(-Rng_comp(i))-H0(1);
    H(i,2) = H(i,2)/(-Rng_comp(i))-H0(2);
    H(i,3) = H(i,3)/(-Rng_comp(i))-H0(3);
    PRange(i) = (Rng_comp(i) - Rng_comp0) + p.wave_l1*( Phase_meas.Nhat(i+1) - Phase_meas.Nhat(1) ); 
    if dual_freq
        PRange(num_meas+i,1) = (Rng_comp(i) - Rng_comp0) + ...
            p.wave_wd*( Phase_meas.Nhat_wd(i+1) - Phase_meas.Nhat_wd(1) ); 
    end
end
residual = Meas - PRange;

HH(:,p.X_POS) = H*p.R_ned2ecef;  %% NOTICE the Rotation matrix, transform from ECEF to TP

if dual_freq
    HH = [HH; HH];
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
