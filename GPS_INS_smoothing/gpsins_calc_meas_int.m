% This function is for calculating double differenced measurements and
% covariance matrices, from the validated phase measurement information
% Phase_meas.prn_phase is the prn list
% Phase_meas.num_ph_sv is total number of satellites with valid phase meas
% Phase_meas.num_ph_meas is the total number of valid phase meas.
function [Meas, R_meas] = gpsins_calc_meas_int(meas_temp, Phase_meas)
global Dual_Freq;
dual_freq = Dual_Freq;
R_meas_temp= zeros(Phase_meas.num_ph_sv-1,1);
if ~dual_freq
    Meas = zeros(Phase_meas.num_ph_sv-1,1);
else
    Meas = zeros(Phase_meas.num_ph_sv-1,2);
end

com_prn = Phase_meas.prn_phase(1); % common sv
for i=2:Phase_meas.num_ph_sv
    prn = Phase_meas.prn_phase(i);
    Meas(i-1,1) = meas_temp.Sat_state(prn).sd_phase_l1 - ...
                          meas_temp.Sat_state(com_prn).sd_phase_l1;
    R_meas_temp(i-1) = meas_temp.Sat_state(prn).R_ph + meas_temp.Sat_state(com_prn).R_ph;
    
    if dual_freq
        Meas(i-1,2) = meas_temp.Sat_state(prn).sd_phase_l2 - ...
                          meas_temp.Sat_state(com_prn).sd_phase_l2;
    end
end

if ~dual_freq
    R_meas = diag(R_meas_temp);
else
    R_meas = diag([R_meas_temp; 5.7^2*R_meas_temp]);
    Meas = [Meas(:,1);Meas(:,2)];
end