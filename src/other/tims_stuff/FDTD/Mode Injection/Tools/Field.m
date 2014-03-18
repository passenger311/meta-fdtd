
function Field(w_in,yst,yfi,xst,dx,type)

% Initialise structure
AndreasHighDielectric = Structure;
AndreasHighDielectric.layers = {            % Define Layers
    'SiO2',...
    'Dielectric(12.25,1)',...
    'SiO2',...
    };

if type == 1
    mode = 'TE';
else
    mode = 'TM';
end

AndreasHighDielectric.d = [125]*1e-9;   % Layer thicknesses
AndreasHighDielectric.pola = mode;          % Light polarization
AndreasHighDielectric.bnd = 0;              % Mode type

% Read Dispersion data
dispfile = 'TM_Mode.txt';
dispdata = dlmread(dispfile);
beta = dispdata(:,1);
realomega = dispdata(:,3);
imagomega = dispdata(:,4);
omega = complex(realomega,imagomega);

% Find nearest omega value
b_in = spline(realomega,beta,w_in);
w_in = spline(beta,omega,b_in);
fprintf(1,sprintf('real omega = %g',real(w_in)))
fprintf(1,sprintf('\n\n imag omega = %g',imag(w_in)))
w = real(w_in);
k0 = w/AndreasHighDielectric.siC;
Neff = b_in/k0;

% Write Neff to file
a = 'Neff.txt';
file_out = fopen(a,'wt');
fprintf(file_out,'%g',Neff);
fclose(file_out);

% Calculate fields
dx = dx*1e-6;
st = (yst+1)*dx;
fi = (yfi-1)*dx + sum(AndreasHighDielectric.d);

fields = AndreasHighDielectric.field_plotter(w_in,b_in,st,-fi,1100);

% Match to FDTD grid
y = (yfi:1:yst)*dx;
ytmp = fields(:,1);

if type == 1
    Etmp = real(fields(:,2));
    Htmp = 299792458*real(fields(:,3));
else
    Htmp = real(fields(:,2));
    Etmp = 299792458*real(fields(:,3));
end

Hf = spline(ytmp,Htmp,y);
Ef = spline(ytmp,Etmp,y);
div = max(abs(Hf));
Hf = Hf./div;
Ef = Ef./div;

% Plot Fields
close all
figure(1)
plot(ytmp,Htmp)
hold on
plot(ytmp,Etmp,'r')
figure(2)
plot(y/dx,Hf)
hold on
plot(y/dx,Ef,'r')
pause(10)

% Write data to file
a = 'source_prism_data.in';
file_output = fopen(a,'wt');

fprintf(file_output,'(SET\n');
fprintf(file_output,'%5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n',xst,xst,1,yfi,yst,1,0,0,1);

for p = 1:length(Ef);
    fprintf(file_output,'%12.6E %12.6E\n',Ef(p), Hf(p));
end;
fprintf(file_output,')SET');
fclose(file_output);

end








