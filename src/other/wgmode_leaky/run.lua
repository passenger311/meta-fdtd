lambdainv=pinv_wl
neffguess=nrefr
fig_on=1
imagpart=nrefr_im;

fh = io.open("run.m","w")
cmd = "[Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( 'geo_injection',"
   ..tostring(lambdainv)..",1,1,"..tostring(neffguess).."+i*"..tostring(imagpart)..","..tostring(nrefr)..","..tostring(fig_on)..
   ","..tostring(cyl_rad)..","..tostring(conv)..");"
fh:write(cmd.."\n")

if silent then
cmd = "quit;"
fh:write(cmd.."\n")
end
fh:close()
os.execute("matlab -r run -nosplash -nodesktop")

