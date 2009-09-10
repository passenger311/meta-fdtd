function saveRCS2 (file, theta, phi, lambda, RCS);

fprintf(1,'\nwriting files: RCS2\n');

max_phi = size(phi,2);
max_theta = size(theta,2);

for n_theta = 1:max_theta
   for n_phi = 1:max_phi 
      filename = sprintf('%sRCS2_%i.dat',file,phi(n_phi)-90);
      fid = fopen(filename, 'a');
      if (fid==-1) error('unable to write file'); end;
      fprintf(fid,'%e\t%e\n',lambda, RCS(n_theta,n_phi));
      fclose(fid);
   end
end
