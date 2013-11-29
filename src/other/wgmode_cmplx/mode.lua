lambdainv=pinv_wl
neffguess=nrefr2
fig_on=1

fh = io.open("run.m","w")
cmd = "[Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( 'geo_injection',"
   ..tostring(lambdainv)..",1,1,"..tostring(neffguess).."+i*0.05,"..tostring(nrefr)..","..tostring(fig_on)..
   ","..tostring(cyl_rad)..");"
fh:write(cmd.."\n")
if silent then
cmd = "quit;"
fh:write(cmd.."\n")
end
fh:close()
os.execute("matlab -r run -nosplash -nodesktop")

