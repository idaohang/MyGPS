
function [xhat, err_flag] = gpsins_correct_xhat_shift(p, xhat, shift)

% position error state
xhat.r_tb_t = xhat.r_tb_t + shift;
