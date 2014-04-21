% This function is to reset the xhat
function [xhat]=gpsins_reset_xhat(start_time, pos)
% needs improvement
cphi = 1;%cos(roll);
sphi = 0;%sin(roll);
ctheta = 1;%cos(pitch);
stheta = 0;%sin(pitch);
cpsi = 1;%cos(heading);
spsi = 0;%sin(heading);

R_t2b = [ ctheta*cpsi                    ctheta*spsi                     -stheta;
          cpsi*stheta*sphi-cphi*spsi     cphi*cpsi+stheta*sphi*spsi      ctheta*sphi;
          cphi*cpsi*stheta+sphi*spsi     -cpsi*sphi+cphi*stheta*spsi     ctheta*cphi ];

xhat.t = start_time;
xhat.R_b2t = R_t2b';
xhat.r_tb_t = pos; %zeros(3,1);
xhat.v_tb_t = zeros(3,1);
xhat.ba_b = zeros(3,1);
xhat.bg_b = zeros(3,1);