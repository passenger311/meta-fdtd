
print("---- geometry")


funcs = require "funcs"
round = funcs.round

infty = 10000000

--- geometrical parameters in real units (um)

real = {}

real.hlength_cpl  = 2                  -- length
real.hwidth_cpl  = 0.1115              -- cpl half width
real.length_in   = 1                  -- input wg length
real.length_out  = 1                  -- output wg length
real.feed_in     = 1.1*real.hwidth_cpl
real.feed_out    = 1.1*real.hwidth_cpl
real.width_in    = 0.08
real.width_out   = real.width_in
real.pad_metal   = 0.02
real.offset_in   = -0.3
real.offset_out   = 0.3
real.dist_inout  = real.offset_out - real.offset_in


real.min_x       = - real.length_in  - real.hwidth_cpl
real.max_x       =   real.length_out + real.hwidth_cpl
real.min_y       = - real.hlength_cpl
real.max_y       =   real.hlength_cpl


--- convert real lengths into grid units

grid = funcs.scale(real, real_dx, true)

---- draw structure

scene_eps = Scene{ value = eps.air }
scene_met = Scene{ value = 0. }

---- draw epsilon structure


cpl = Box{ from={ -grid.hwidth_cpl, -grid.hlength_cpl-infty }, to={ grid.hwidth_cpl, grid.hlength_cpl+infty } }

oy = grid.offset_in - grid.width_in/2
ox = - grid.hwidth_cpl + grid.feed_in

cin = Box{ from={ ox, oy }, to ={ ox-infty, oy+grid.width_in } }

oy = grid.offset_out - grid.width_out/2
ox = grid.hwidth_cpl - grid.feed_in

cout = Box{ from={ ox, oy }, to ={ ox+infty, oy+grid.width_out } }

scene_eps:add{ cpl, value = eps.si }
scene_eps:add{ cin, value = eps.sio2 }
scene_eps:add{ cout, value = eps.sio2 }


---- draw metal structure

m_cpl = Box{ from={ -grid.hwidth_cpl-grid.pad_metal, -grid.hlength_cpl-infty }, to={ grid.hwidth_cpl+grid.pad_metal, grid.hlength_cpl+infty } }

oy = grid.offset_in - grid.width_in/2
ox = - grid.hwidth_cpl 

m_cin = Box{ from={ ox, oy-grid.pad_metal }, to ={ ox-grid.length_in-infty, oy+grid.width_in+grid.pad_metal } }

oy = grid.offset_out - grid.width_out/2
ox = grid.hwidth_cpl

m_cout = Box{ from={ ox, oy-grid.pad_metal }, to ={ ox+grid.length_out+infty, oy+grid.width_out+grid.pad_metal } }


scene_eps:add{ BinaryAndNot{BinaryAndNot{BinaryAndNot{m_cpl, cpl},cin},cout} , value = -1e20 }
scene_eps:add{ BinaryAndNot{m_cin, cin} , value = -1e20 }
scene_eps:add{ BinaryAndNot{m_cout, cout} , value = -1e20 }


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


--- end of geo

