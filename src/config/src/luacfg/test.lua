
sc1 = geo.Scene();
sc2 = geo.Scene();
x = geo.Sphere{pos={1.,2.,3.},rad=4.7}


for i=1,2000,1 do
   x =  geo.Sphere{at={1.,2.,3.},rad=4.7}
   y =  geo.Sphere{at={1.,2.,3.},rad=4.7}
   z =  geo.BinaryOr{x,y}
   sc1:add(x)
   sc1:add(y)
   sc1:add(z)
end

for i=1,2000000,1 do
   x =  geo.Sphere{at={1.,2.,3.},rad=4.7}
   r = geo.Transform{x, move={7,3,2}}
   sc2:add(r)
end