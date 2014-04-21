function [num_meas, Meas, R_meas, prn_dd] = gpsins_calc_meas_msc(meas_temp)

num_meas = 0;
Meas_temp = zeros(meas_temp.num_sats,1);
R_meas_temp = zeros(meas_temp.num_sats,1);
prn_dd_temp = zeros(meas_temp.num_sats,2);

if meas_temp.num_dgps>0
    for i=1:meas_temp.num_sats
        prn = meas_temp.prnlist(i);
        if meas_temp.Sat_state(prn).sat_valid
            if (meas_temp.Sat_state(prn).dgps_age >=0 && meas_temp.Sat_state(prn).dgps_age <=5)
                com_prn_ph = prn;      % set ph common sat
                break
            end
        end
    end
end

for i=1:meas_temp.num_sats
    prn = meas_temp.prnlist(i);
    if meas_temp.Sat_state(prn).sat_valid
        if prn~=com_prn_ph && (meas_temp.Sat_state(prn).dgps_age >=0 && meas_temp.Sat_state(prn).dgps_age <=5)
            num_meas = num_meas+1;
            prn_dd_temp(num_meas,:) = [prn, com_prn_ph];
            Meas_temp(num_meas) = (meas_temp.Sat_state(prn).sd_phase_l1 - meas_temp.Sat_state(com_prn_ph).sd_phase_l1);
            R_meas_temp(num_meas) = meas_temp.Sat_state(prn).R_ph + meas_temp.Sat_state(com_prn_ph).R_ph;
        else
            disp(['At time ', num2str(meas_temp.imu_tm),' for Sv ', num2str(prn), ' no measurement used']);
        end
    else
        disp(['At time ', num2str(meas_temp.imu_tm),' for Sv ', num2str(prn), ' no measurement used. Since invalid.']);
    end
end

R_meas = diag(R_meas_temp(1:num_meas));
prn_dd = prn_dd_temp(1:num_meas,:);
Meas = Meas_temp(1:num_meas);