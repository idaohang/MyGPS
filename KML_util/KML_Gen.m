%% Load data
ins20140123_1558;

%% Generate KML
start_time = 3;
end_time = 800;
filename = '1558.kml';

% trim data
i = 1;
while ( ins(i,2)<start_time )
    i=i+1;
end
j=i+1;
while ( ins(j,2)<end_time )
    j=j+1;
    if j==length(ins(:,1))
        break;
    end
end
pos_tp = ins(i:j,3:5);
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

pos_lon_lat = [pos_lla(:,2),pos_lla(:,1)];

myKML_write(filename,pos_lla);

%kmlwrite('1558.kml',pos_lla(:,1),pos_lla(:,2),pos_lla(:,3));
