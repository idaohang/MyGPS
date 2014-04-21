%% This function is to use Doppler measurement, and also decouple the noise
function [HH_2, HH, R_meas, residual] = gps_doppler_output(p, rpos_ecef, rpos_ecef_2, meas_temp, meas_temp_2)
%%
[num_meas, Meas, R_meas, prn_dd] = gpsins_calc_meas(meas_temp);
[num_meas_2, Meas_2, R_meas_2, prn_dd_2] = gpsins_calc_meas(meas_temp_2);

svpos_ecef = zeros(3,1);
H0 = zeros(1,3);
H_2 = zeros(num_meas,3);
H = zeros(num_meas,3);
Rng_comp_2 = zeros(num_meas,1);
Rng_comp = zeros(num_meas,1);
PRange_2 = zeros(num_meas,1);
PRange = zeros(num_meas,1);
HH = zeros(num_meas,p.X_STATES);
HH_2 = zeros(num_meas,p.X_STATES);

for i=1:num_meas
    prn0 = prn_dd(i,2);
    svpos_ecef = meas_temp.Sat_state(prn0).sv_pos_ecef;
    H0(1) = svpos_ecef(1) - rpos_ecef(1);
    H0(2) = svpos_ecef(2) - rpos_ecef(2);
    H0(3) = svpos_ecef(3) - rpos_ecef(3);
    Rng_comp0 = norm( H0 );
    H0(1) = H0(1)/(-Rng_comp0);
    H0(2) = H0(2)/(-Rng_comp0);
    H0(3) = H0(3)/(-Rng_comp0);
    
    prni = prn_dd(i,1);
    svpos_ecef = meas_temp.Sat_state(prni).sv_pos_ecef;
    H(i,1) = svpos_ecef(1) - rpos_ecef(1);
    H(i,2) = svpos_ecef(2) - rpos_ecef(2);
    H(i,3) = svpos_ecef(3) - rpos_ecef(3);
    Rng_comp(i) = norm( H(i,:) );
    H(i,1) = H(i,1)/(-Rng_comp(i))-H0(1);
    H(i,2) = H(i,2)/(-Rng_comp(i))-H0(2);
    H(i,3) = H(i,3)/(-Rng_comp(i))-H0(3);
    PRange(i) = (Rng_comp(i) - Rng_comp0);
    
    %%
    svpos_ecef_2 = meas_temp_2.Sat_state(prn0).sv_pos_ecef;
    H0(1) = svpos_ecef_2(1) - rpos_ecef_2(1);
    H0(2) = svpos_ecef_2(2) - rpos_ecef_2(2);
    H0(3) = svpos_ecef_2(3) - rpos_ecef_2(3);
    Rng_comp0_2 = norm( H0 );
    H0(1) = H0(1)/(-Rng_comp0);
    H0(2) = H0(2)/(-Rng_comp0);
    H0(3) = H0(3)/(-Rng_comp0);
    
    prni = prn_dd(i,1);
    svpos_ecef_2 = meas_temp_2.Sat_state(prni).sv_pos_ecef;    
    H_2(i,1) = svpos_ecef_2(1) - rpos_ecef_2(1);
    H_2(i,2) = svpos_ecef_2(2) - rpos_ecef_2(2);
    H_2(i,3) = svpos_ecef_2(3) - rpos_ecef_2(3);
    Rng_comp_2(i) = norm( H_2(i,:) );
    H_2(i,1) = H_2(i,1)/(-Rng_comp_2(i))-H0(1);
    H_2(i,2) = H_2(i,2)/(-Rng_comp_2(i))-H0(1);
    H_2(i,3) = H_2(i,3)/(-Rng_comp_2(i))-H0(1);
    PRange_2(i) = (Rng_comp_2(i) - Rng_comp0_2);
end
residual = Meas - Meas_2 - PRange + PRange_2;
HH(:,p.X_POS) = H*p.R_ned2ecef;  %% NOTICE the Rotation matrix, transform from ECEF to TP
HH_2(:,p.X_POS) = H_2*p.R_ned2ecef;  %% NOTICE the Rotation matrix, transform from ECEF to TP

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
