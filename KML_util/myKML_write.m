function flag = myKML_write(filename,pos_lla)
fid=fopen(filename, 'w');
header = '<?xml version="1.0" encoding="utf-8"?>\n';
fprintf(fid, header);
fprintf(fid, '<kml> xmlns="http://www.opengis.net/kml/2.2"\n');
fprintf(fid, '  <Document>\n');
fprintf(fid, '    <Style id="DGPS">\n');
fprintf(fid, '      <LineStyle>\n');
fprintf(fid, '        <color>ff0000ff</color>\n');
fprintf(fid, '        <width>4</width>\n');
fprintf(fid, '      </LineStyle>\n');
fprintf(fid, '      <PolyStyle>\n');
fprintf(fid, '        <color>ff0000ff</color>\n');
fprintf(fid, '      </PolyStyle>\n');
fprintf(fid, '    </Style>\n');
fprintf(fid, '    <Placemark>\n');
fprintf(fid, '      <Style>\n');
fprintf(fid, '        <IconStyle>\n');
fprintf(fid, '          <Icon>\n');
fprintf(fid, '            <href>root://icons/palette-3.png</href>');
fprintf(fid, '            <x>0</x>\n');
fprintf(fid, '            <y>0</y>\n');
fprintf(fid, '            <w>32</w>\n');
fprintf(fid, '            <h>32</h>\n');
fprintf(fid, '          </Icon>\n');
fprintf(fid, '        </IconStyle>\n');
fprintf(fid, '      </Style>\n');
fprintf(fid, '      <name>Marker 1</name>\n');
fprintf(fid, '      <description>UCR CERT, Center for Environmental Research and Technology</description>\n');
fprintf(fid, '      <Point>\n');
coordinates0 = ['      <coordinates>',num2str(pos_lla(1,2),16),',',num2str(pos_lla(1,1),16),',0</coordinates>\n'];
%fprintf(fid, '      <coordinates>-117.335583209000,34.0005672342000,0</coordinates>\n');
fprintf(fid,coordinates0);
fprintf(fid, '      </Point>\n');
fprintf(fid, '    </Placemark>\n');

fprintf(fid, '    <Placemark>\n');
fprintf(fid, '      <LookAt>\n');
longitude0 = ['        <longitude>', num2str(mean(pos_lla(:,2)),16),'</longitude>\n'];
latitude0  = ['        <latitude>', num2str(mean(pos_lla(:,1)),16),'</latitude>\n'];
altitude0  = ['        <altitude>', num2str(mean(pos_lla(:,3)),16),'</altitude>\n'];
fprintf(fid, longitude0);
fprintf(fid, latitude0);
fprintf(fid, altitude0);
fprintf(fid, '      </LookAt>\n');
fprintf(fid, '      <styleUrl>#DGPS</styleUrl>\n'); 
fprintf(fid, '      <LineString>\n');
fprintf(fid, '        <coordinates>\n');
for i = 1:length(pos_lla(:,1))
    coordinate = ['          ',num2str(pos_lla(i,2),16),',',num2str(pos_lla(i,1),16),'\n'];
    fprintf(fid, coordinate);
end
fprintf(fid, '        </coordinates>\n');
fprintf(fid, '      </LineString>\n');
fprintf(fid, '    </Placemark>\n');
fprintf(fid, '  </Document>\n');
fprintf(fid, '</kml>\n');




