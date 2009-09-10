function mysurf2(number)

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

d = {'cap','capnp'};
for l=1:2
  dirname = sprintf('task_%s%i/',d{l},number)
  cd(dirname);
  for j = 1:km
    filename = sprintf('RCS_%i.dat',j);
    fid1=fopen(filename,'r');
    k = 1;
    tline = fgetl(fid1);
    while feof(fid1) == 0
      tline = fgetl(fid1);
      tmp = sscanf(tline, '%e %e');
      RCS(j,k,l,:) = tmp; 
      k = k+1;
    end;
    fclose(fid1);
  end;
  cd ..;
end;

angle = linspace(-180,180,361);
figure(1)
surf(angle,lambda(20:101),RCS(20:101,:,1,2))
shading interp
figure(2)
surf(angle,lambda(20:101),RCS(20:101,:,2,2))
shading interp
figure(3)
surf(angle,lambda(20:101),RCS(20:101,:,2,2)-RCS(20:101,:,1,2))
shading interp
