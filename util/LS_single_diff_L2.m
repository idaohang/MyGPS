function [rpos_ecef, delta_t] = LS_single_diff_L2(meas_temp, init_pos_ecef)
%% options
DGPS_ONLY = 0;
CODE_ONLY = 0;
%% constants
NUM_SATS = 15;
NUM_PRNS = 32;
GPS_TH = 0.0001;
%%
rpos_ecef = zeros(4,1);
rpos_ecef(1:3)=init_pos_ecef;
svpos_ecef = zeros(4,1);
delta = zeros(4,1);
H = zeros(meas_temp.num_sats,4);
Rng_comp = zeros(meas_temp.num_sats,1);
PRange = zeros(meas_temp.num_sats,1);
Meas = zeros(meas_temp.num_sats,1);
R_meas = zeros(meas_temp.num_sats,1);
for i=1:meas_temp.num_sats
    prn = meas_temp.prnlist(i);
    if meas_temp.Sat_state(prn).sat_valid
        if meas_temp.Sat_state(prn).locked && (~DGPS_ONLY)
            Meas(i) = meas_temp.Sat_state(prn).ph_l2_rng;
            R_meas(i) = meas_temp.Sat_state(prn).R_ph;
        else if ( meas_temp.Sat_state(prn).dgps_age >=0 && meas_temp.Sat_state(prn).dgps_age <=5) && (~CODE_ONLY)
                Meas(i) = meas_temp.Sat_state(prn).sd_code_l2;
                R_meas(i) = meas_temp.Sat_state(prn).R_cd;
                disp(['At time ', num2str(meas_temp.imu_tm),' Sv ', num2str(prn), ' use sd_code_l2']);
            else
                Meas(i) = meas_temp.Sat_state(prn).code_l2;
                R_meas(i) = meas_temp.Sat_state(prn).R_cd;
                disp(['At time ', num2str(meas_temp.imu_tm),' Sv ', num2str(prn), ' use code_l2']);
            end
        end
    else
        % error print said this sat at this time is not valid
    end
end

recur_count = 0;
while(recur_count==0 || norm(delta)> 0.01*GPS_TH)
    for i=1:meas_temp.num_sats
        prn = meas_temp.prnlist(i);
        svpos_ecef = meas_temp.Sat_state(prn).sv_pos_ecef;
        H(i,1) = svpos_ecef(1) - rpos_ecef(1);
        H(i,2) = svpos_ecef(2) - rpos_ecef(2);
        H(i,3) = svpos_ecef(3) - rpos_ecef(3);
        Rng_comp(i) = norm( H(i,1:3) );
        H(i,1) = H(i,1)/(-Rng_comp(i));
        H(i,2) = H(i,2)/(-Rng_comp(i));
        H(i,3) = H(i,3)/(-Rng_comp(i));
        H(i,4) = 1;
        PRange(i) = (Rng_comp(i) + rpos_ecef(4));
    end
    Cov_meas = diag(1./R_meas);
    residual = Meas - PRange;
    delta = (H'*Cov_meas*H)\(H'*Cov_meas*residual);  % inv(H'*inv(Cov_meas)*H)*H'*Cov_meas*residual
    rpos_ecef = rpos_ecef + delta;
    recur_count=recur_count+1;
end

delta_t = rpos_ecef(4);
rpos_ecef = rpos_ecef(1:3);

end
