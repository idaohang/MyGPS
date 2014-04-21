function p_ned = convert_ecef2NED(p_ecef)

%base ecef
ecef_p_b(1) = -2691969.038;
ecef_p_b(2) = -4282611.458;
ecef_p_b(3) = 3871940.615;

geod = ecef2geod(ecef_p_b);

% store lat, lon, alt (in deg)
lat = dms2deg(geod(1,:));
lon = dms2deg(geod(2,:));
alt = geod(3,1);


% all in rad
% lontitude: lambda =
% latitude: phi =
lambda = -lon*pi/180;
phi = lat*pi/180;
R_ecef2T = [-sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi);
    -sin(lambda) cos(lambda) 0;
    -cos(phi)*cos(lambda) -cos(phi)*sin(lambda) -sin(phi)];


p_ned = (R_ecef2T * (p_ecef' - ecef_p_b'))';


function [y]=dms2deg(x)
if x(1)>0,
    y = x(1)+(x(2) + x(3)/60)/60;
else
    y = x(1)-(x(2) - x(3)/60)/60;
end