
-- OPENMP benchmark script
-- J. Hamm 2013

local parallel = 0
local tot_time1
for threads=1,2 do 
   -- run config file
   local cmd = "./luacfg ./run_config.lua "..tostring(100).." 48000"
   os.execute(cmd.." >luacfg.log")
   -- run meta
   cmd = "OMP_NUM_THREADS="..tostring(threads).." ./meta"
   os.execute(cmd.." >meta.log")
   -- get total fdtd loop time from log
   file = assert(io.open("meta.log"))
   output = file:read('*all')
   file:close()
   tot_time = tonumber(string.match(output,"%(TIMER%)%s+Total%s+:%s+(.*)%s+secs"))
   -- calculate speedup
   if not tot_time1 then 
      tot_time1 = tot_time 
   else 
      parallel = threads * ( 1 - tot_time / tot_time1 ) / ( threads - 1 ) 
   end 
   local speedup = tot_time1/tot_time
   -- calculate parallel percentage
   
   -- print timing
   print(threads, tot_time, speedup, parallel )
end
