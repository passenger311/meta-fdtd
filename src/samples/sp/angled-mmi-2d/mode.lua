
lambdainv = select(1,...)
neffguess = select(2,...)
silent = select(3,...)

fh = io.open("run.m","w")
cmd = "[ Ex,Ey,Hx,Hy,range] = wgmode3d_2( 'geo_inj.in',"
   ..tostring(lambdainv)..",1,1,"..tostring(neffguess)..");"
fh:write(cmd.."\n")
cmd = "savetfsfinj_2( 'tfsfex.set', Ex, Hy, range, range );"
fh:write(cmd.."\n")
cmd = "savetfsfinj_2( 'tfsfey.set', Ey, Hx, range, range );"
fh:write(cmd.."\n")
if silent then
cmd = "quit;"
fh:write(cmd.."\n")
end
fh:close()

os.execute("matlab -r run -nosplash -nodesktop")
