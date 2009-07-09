module ("funcs", package.seeall)

infty = 100000

function nmax(eps)
   local epsmax = 0
   for k,v in pairs(eps) do
      epsmax = math.max(epsmax,v)
   end
   return math.sqrt(epsmax)
end


function round(f)
   return math.floor(f+0.5)
end


function scale(tab, real_dx, talk)
   local stab = {}
   local function _scale(tab, ntab, real_dx, str)
      local kstr
      for k,v in pairs(tab) do
	 if type(k) == 'string' then 
	    kstr = "['"..tostring(k).."']" 
	 else
	    kstr = "["..tostring(k).."]" 
	 end
	 if type(v) == 'table' then 
	    ntab[k] = {}
	    _scale(v, ntab[k], real_dx, str..kstr)
	 elseif type(v) == 'number' then
	    if talk then
	       print("scale->"..str..kstr, v, v/real_dx)
	    end
	    ntab[k] = v/real_dx
	 end
      end
   end
   _scale(tab,stab,real_dx, "")
   return stab
end


function offset(tab)
   local off = {}
   local l = 0
   for i,v in ipairs(tab) do
      off[i] = l
      l = l + v
   end
   off[0] = off[1]-infty
   off[#tab+1] = l 
   off[#tab+2] = l+infty
   off.total = l
   return off
end   
