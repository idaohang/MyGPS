% This script plots rover trajectory with satellite number (in) color 
% on Google Map. The function plot_google_map() is called to use GMap API

sat_num_log20140130_1134;

pos_tp = satnumlog(:,3:5);
sat_num = satnumlog(:,6);
buff_size = length( pos_tp(:,1) );
pos_ecef = zeros(buff_size, 3);

p.base_ecef = [-2430691.399 -4704193.471 3544329.054]';
p.base_lla = ecef2lla(p.base_ecef');  % WGS84 datum

% p.R_ned2ecef = [0.2565316863   0.8884106506  0.3806809816
%               0.4964724789  -0.4590495789  0.7367418556
%               0.8292807556   0            -0.5588322005];

lambda = p.base_lla(2)*pi/180;
phi = p.base_lla(1)*pi/180;
p.R_ecef2ned = [-sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi);
                    -sin(lambda)          cos(lambda)          0;
                -cos(phi)*cos(lambda) -cos(phi)*sin(lambda) -sin(phi)];

for i=1:buff_size
    pos_ecef(i,:) = (p.R_ecef2ned'*pos_tp(i,:)')+ p.base_ecef;
end

pos_lla = ecef2lla(pos_ecef);  % WGS84 datum

lat = pos_lla(:,1);
lon = pos_lla(:,2);

figure(30)
clf;
scatter(lon, lat, 50, sat_num, 'filled'); colorbar;
plot_google_map()

figure(31)
clf;
scatter(lon, lat, 50, sat_num, 'filled'); colorbar;
plot_google_map('MapType', 'satellite')

figure(32)
clf;
scatter(lon, lat, 50, sat_num, 'filled'); colorbar;
plot_google_map('MapType', 'hybrid')