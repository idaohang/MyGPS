function  Phase_meas = validate_phase_meas(nav, data, Xhat)
global Dual_Freq;
dual_freq = Dual_Freq;
p = nav.param;
K = length(data);

% FIXME: improve this
prn_phase = merge_prn_list(data);
%prn_phase = [13;7;6;10;19;8;20;23;16];
[num_ph_sv, ~] = size(prn_phase);
[num_epoch, ~] = size(Xhat);
num_ph_meas = (num_ph_sv-1)*( num_epoch-1 ); % double difference
Nhat = zeros(num_ph_sv,1);
if dual_freq
    Nhat_wd = zeros(num_ph_sv,1);
end

for i = 1:num_ph_sv
    prn = prn_phase(i);
    for j = 1:size(Xhat)-1
        sv_pos = data.gps(j,1).Sat_state(1,prn).sv_pos_ecef'; 
        r_pos= p.ecef_p_b + p.R_ned2ecef*Xhat(j+1,1).r_tb_t;
        sd_phase = data.gps(j,1).Sat_state(1,prn).sd_phase_l1;
        Nhat(i) = Nhat(i) + sd_phase*p.wave_l1 - norm( r_pos - sv_pos );
        if dual_freq
            sd_phase_wd = data.gps(j,1).Sat_state(1,prn).sd_phase_l1...
                          - data.gps(j,1).Sat_state(1,prn).sd_phase_l2;
            Nhat_wd(i) = Nhat_wd(i) + sd_phase_wd*p.wave_wd - norm( r_pos - sv_pos );
        end
    end
    Nhat(i) = round( Nhat(i)/p.wave_l1/j );
    if dual_freq
        Nhat_wd(i) = round( Nhat_wd(i)/p.wave_wd/j );
    end
end

Phase_meas.prn_phase = prn_phase;
Phase_meas.num_ph_sv = num_ph_sv;
Phase_meas.num_ph_meas = num_ph_meas;
Phase_meas.Nhat = Nhat;
if dual_freq
    Phase_meas.Nhat_wd = Nhat_wd;
end