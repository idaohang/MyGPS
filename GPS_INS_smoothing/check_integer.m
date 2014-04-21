function [N] = check_integer(Phase_meas,gps)
N = zeros(Phase_meas.num_ph_sv,1);
for i = 1:Phase_meas.num_ph_sv
    prn = Phase_meas.prn_phase(i);
    N(i) = gps.Sat_state(prn).N_l1;
end