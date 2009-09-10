function mysurf(dir,number)
% mysurf(dir,number)
%
% The function mysurf reads the differential scattering cross section for different wavelengths from RCS_%i.dat and compiles these in one graph
% 
% Input:
% - dir = string of directory name %s: task_%s%i, i.e. task_cap1
% - number = %i in above directory name
%

fid1 = fopen('invlambda.in','r');                            
fid2 = fopen('data.save','r');                               
k = 1;                                                       
while feof(fid1) == 0                                        
  tline = fgetl(fid1);                                          
  invlambda=sscanf(tline,'%e');                                 
  lambda(k) = 1/invlambda;                                      
  k = k+1;                                                      
end;                                                            
km = k-1;                                                    
tline = fgetl(fid2);                                         
dx = sscanf(tline,'%e');                                     
                                                             
                                                             
fclose(fid1);                                                
fclose(fid2);

dirname = sprintf('task_%s%i/',dir,number)
cd(dirname);
for j = 1:km
  filename = sprintf('RCS_%i.dat',j);
  fid1=fopen(filename,'r');
  k = 1;
  tline = fgetl(fid1);
  while feof(fid1) == 0
    tline = fgetl(fid1);
    tmp = sscanf(tline, '%e %e');
    RCS(j,k,:) = tmp; 
    k = k+1;
  end;
  fclose(fid1);
end;
cd ..;

angle = linspace(-180,180,361);
surf(angle,lambda(20:101),RCS(20:101,:,2))
shading interp

