
sc1 = geo.Scene();
sc2 = geo.Scene();
x = geo.Sphere{pos={1.,2.,3.},rad=4.7}


for i=1,2000,1 do
   x =  geo.Sphere{at={1.,2.,3.},radius=4.7}
   y =  geo.Box{at={1.,2.,3.},size={10,10,10}}
   z =  geo.BinaryOr{x,y}
   sc1:add(x)
   sc1:add(y)
   sc1:add(z)
end

for i=1,2000000,1 do
   x =  geo.Sphere{at={1.,2.,3.},radius=4.7}
   y =  geo.Ellipsoid{at={1.,2.,3.},size={7,9,10}}
   y =  geo.Cylinder{at={1.,2.,3.},radius=8,height=4}
   r = geo.Transform{x, move={7,3,2}}
   sc2:add(r)
end