clear all

fidsave=fopen('backscat.dat','w');
fidsave3=fopen('backscatcrosssec.dat','w');
fidsave2=fopen('forwardscat.dat','w');
fidsave4=fopen('forwardscatcrosssec.dat','w');
for i = 500:999
   filename=sprintf('%i.out',i);
   fprintf(1,'loading file: %s\n', filename);
   fid = fopen(filename, 'r');
   if (fid==-1) error('file does not exist'); end;
   tmp = textscan(fid, '%n %n %n %n %n %n %n');
   fclose(fid);
   angle = tmp{1}; sizemax= size(angle,1);
   sumparaperp = zeros(sizemax,1);para = zeros(sizemax,1);
   perp = zeros(sizemax,1);

   sumparaperp = tmp{2}/(2*pi/(i*1d-9))^2;
   para = tmp{7}/(2*pi/(i*1d-9))^2;
   perp = tmp{4}/(2*pi/(i*1d-9))^2;

   fprintf(fidsave,'%i      %e\n', i,sqrt(sumparaperp(sizemax)/2));
   fprintf(fidsave2,'%i      %e\n', i,sqrt(sumparaperp(1)/2));
   fprintf(fidsave3,'%i      %e\n', i,sumparaperp(sizemax)*(2*pi/(i*1d-9))^2/pi/2);
   fprintf(fidsave4,'%i      %e\n', i,sqrt(sumparaperp(1)*(2*pi/(i*1d-9))^2/pi/2));
end
fclose(fidsave);
fclose(fidsave2);
fclose(fidsave3);
fclose(fidsave4);
