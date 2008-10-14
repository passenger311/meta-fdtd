
x = geo.Sphere{pos={1.,2.,3.},rad=4.7}


for i=1,2000,1 do
   x =  geo.Sphere{at={1.,2.,3.},rad=4.7}
   y =  geo.Sphere{at={1.,2.,3.},rad=4.7}
   z =  geo.BinaryOr{x,y}
end