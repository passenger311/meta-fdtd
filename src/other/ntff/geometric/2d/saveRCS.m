function saveRCS (file, freqnumber, theta, phi, RCS);

filename = sprintf('%sRCS_%i.dat',file,freqnumber);
fprintf(1,'\nwriting file: %s\n', filename);
fid = fopen(filename, 'w');
if (fid==-1) error('unable to write file'); end;
max_phi = size(phi,2);
max_theta = size(theta,2);

n_theta=1;
for n_phi=1:max_phi
  fprintf(fid, '\n%e', phi(n_phi)-90);
%  for n_phi=1:max_phi
    fprintf(fid,'       %e',RCS(n_theta,n_phi));
%  end
end
fprintf(fid,'\n');
fclose(fid);
