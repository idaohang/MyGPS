function [H, R_meas, residual] = gpsins_gps_output_loose(p, rpos_ecef, meas_temp, dgps)
%% options
DGPS_ONLY = dgps;

[pos_ecef_gps] = LS_double_diff(meas_temp, rpos_ecef, DGPS_ONLY, 1);

H = zeros(3,p.X_STATES);
H(:,p.X_POS) = eye(3)*p.R_ned2ecef;
residual = pos_ecef_gps - rpos_ecef;
if dgps
    R_meas = diag([0.1,0.1,0.5].^2);
else
    R_meas = diag([0.01,0.01,0.05].^2);
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
