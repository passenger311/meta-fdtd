
print("---- geometry")


funcs = require "funcs"
round = funcs.round

infty = 10000000

real = {}

strip = {}
strip_eps = {}
real.strip = strip

function addstrip( tab )
   table.insert(strip, tab)
   table.insert(strip_eps,tab.eps)
end


--- geometrical parameters in real units (um)

D = 0.0
 
real.length_c   = 1.60
real.length_1   = 0.50
real.length_3   = 0.15
real.length_2   = 0.05

real.min_x      = -0.5
real.max_x      = real.length_c + 0.5
real.min_y      = 0
real.max_y      = 2.
real.width_wg    = 0.25
real.pad_m       = 0.2

w = real.width_wg  
LC = real.length_c
L2 = real.length_2
L1 = real.length_1
L3 = real.length_3

--- strip definition

addstrip{
   at =  { 0, -infty },
   to_y = 1 + L1,
   eps = eps.sio2
}

addstrip{
   at =  { LC, 1 - L3 },
   to_y = infty,
   eps = eps.sio2
}

addstrip{
   at =  { -L2, 1.0 },
   to_x = -w*D,
   eps = eps.si
}

addstrip{
   at =  { w*D , 1.0 },
   to_x = LC-w*D,
   eps = eps.si
}

addstrip{
   at =  { LC+w*D, 1.0 },
   to_x = LC+L2,
   eps = eps.si
}


--- convert real lengths into grid units

grid = funcs.scale(real, real_dx, true)


--- calculate strip positions

w = grid.width_wg/2  

for i,s in ipairs(grid.strip) do
   if s.to_x then
      s.from = { s.at[1]-w, s.at[2]-w }
      s.to = { s.to_x+w, s.at[2]+w }
   elseif s.to_y then
      s.from = { s.at[1]-w, s.at[2]-w }
      s.to = { s.at[1]+w, s.to_y+w }
   end
end

---- draw structure

scene_eps = Scene{ value = eps.air }

---- draw epsilon structure

p = grid.pad_m

f = Box{from={0,0}, to={0,0} } 

for i,s in ipairs(grid.strip) do
   
   s.box1 = Box{ from=s.from, to=s.to }
   s.box2 = Box{ from={s.from[1]-p, s.from[2]-p}, to={ s.to[1]+p, s.to[2]+p }}
   s.frame = BinaryAndNot{ s.box2, s.box1 }

   scene_eps:add{ s.box1, value = strip_eps[i] }

   f = BinaryOr{ f, s.box1 }

   scene_eps:add{ BinaryAndNot{ s.frame,f }, value = eps.ag }


end


---- calculate grid coordinates

imin = round( grid.min_x )
imax = round( grid.max_x )
jmin = round( grid.min_y )
jmax = round( grid.max_y )

kmin = 0
kmax = 0

pad1=0
cpml = 11

imin0 = imin - pad1 - cpml
imax0 = imax + pad1 + cpml
jmin0 = jmin - pad1 - cpml
jmax0 = jmax + pad1 + cpml
kmin0 = kmin
kmax0 = kmax

print("irange (grid) = ", imin, imax)
print("jrange (grid) = ", jmin, jmax)
print("krange (grid) = ", kmin, kmax)

print("irange (real) = ", imin*real_dx, imax*real_dx)
print("jrange (real) = ", jmin*real_dx, jmax*real_dx)
print("krange (real) = ", kmin*real_dx, kmax*real_dx)


--- source definition


src_strip = 1
src_j = round(grid.min_y) + 30

strip = grid.strip[src_strip]
src_i0 = round( strip.at[1]-w)
src_i1 = round( strip.at[1]+w)

src_pad = 20

betaeff = math.sqrt(strip_eps[src_strip])



--- end of geo

