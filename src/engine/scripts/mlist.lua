#!../config/lua/src/lua

arg1 = select(1,...)

if arg1 == "-s" then
   sep = " "
   argc = 2
else
   sep = ","
   argc = 1
end
out = ""

for i = argc, select('#',...) do
   fh = io.open(select(i,...),"r")
   if fh then
      while true do
	 line = fh:read("*line");
	 if line == nil then break end
	 name = line:gsub("%s","")
	 eh = io.open(name..".f90","r")
	 if eh then --- exists?
	    eh:close()
	    if sep == "," then
	       out = out..string.upper(name)..sep
	    else
	       out = out..string.lower(name)..sep
	    end
	 end
      end
	 fh:close()
   end
end

out = out:gsub(",$",""):gsub("%s$",""):gsub("\n$","")
io.write(out)