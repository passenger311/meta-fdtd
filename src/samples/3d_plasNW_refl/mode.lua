lambdainv=pinv_wl
neffguess=nrefr
fig_on=0

fh = io.open("run.m","w")
cmd = "[Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( 'geo_injection',"
   ..tostring(lambdainv)..",1,1,"..tostring(neffguess)..","..tostring(fig_on)..");"
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfex_re.in', real(Ex), real(Hy), range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfex_im.in', imag(Ex), imag(Hy), range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfey_re.in', real(Ey), -real(Hx), range, range );" -- use -Hx because H vector of tfsfinj points in neg x direction
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfey_im.in', imag(Ey), -imag(Hx), range, range );"
fh:write(cmd.."\n")

cmd = "wgsave3k3( sprintf('Emode%2.2f_re.set',neff), real(Ex), real(Ey), real(Ez), range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k3( sprintf('Emode%2.2f_im.set',neff), imag(Ex), imag(Ey), imag(Ez), range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k3( sprintf('Hmode%2.2f_re.set',neff), real(Hx), real(Hy), real(Hz), range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k3( sprintf('Hmode%2.2f_im.set',neff), imag(Hx), imag(Hy), imag(Hz), range, range );"
fh:write(cmd.."\n")

if silent then
cmd = "quit;"
fh:write(cmd.."\n")
end
fh:close()
os.execute("matlab -r run -nosplash -nodesktop")

