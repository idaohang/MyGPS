% This function transform the trajectory X from the data struct to matrix.
function [X_mat] = struct_to_matrix(X, p)
n = length(X);
X_mat = zeros(15,n);

for i = 1:n 
    [roll, pitch, yaw] = dcm2angle( p.R_ned2ecef*X(i,1).R_b2t', 'XYZ' ); % 
    pos_ecef = p.R_ned2ecef*X(i,1).r_tb_t + p.ecef_p_b;
    vel_ecef = p.R_ned2ecef*X(i,1).v_tb_t;
    
    X_mat(:,i) = [ pos_ecef', vel_ecef', roll, pitch, yaw, X(i,1).ba_b', X(i,1).bg_b']';
end
    