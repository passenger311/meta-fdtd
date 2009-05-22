

filename = select(1,...)
lambdainv = select(2,...)
neffguess = select(3,...)
silent = select(4,...)

fh = io.open("run.m","w")
cmd = "[ Ex,Ey,Hx,Hy,range] = wgmode3k( '"..filename.."' ,"
   ..tostring(lambdainv)..",1,1,"..tostring(neffguess)..");"
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfex.set', Ex, Hy, range, range );"
fh:write(cmd.."\n")
cmd = "wgsave3k( 'tfsfey.set', Ey, Hx, range, range );"
fh:write(cmd.."\n")
if silent then
cmd = "quit;"
fh:write(cmd.."\n")
end
fh:close()

os.execute("matlab -r run -nosplash -nodesktop")
