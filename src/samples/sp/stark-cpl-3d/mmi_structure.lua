
taskid = select(1,...)
--- input

l1 = length_wg1
l2 = length_mmi
l3 = length_wg2
w1 = hwidth_mmi
w2 = hwidth_sep
w3 = hwidth_wg

alpha = math.pi/180. * angle

padx = 20
pady = 20

----

y1 = 0 - 100
y2 = l1
y3 = y2+l2
y4 = y3+l3 + 100

x1_1 = w3
x2_1 = x1_1
x2_2 = w1
x3_1 = w2
x3_2 = 2*w3/math.cos(alpha) + x3_1
x3_3 = w1
x4_1 = x3_1 + math.tan(alpha)*l3
x4_2 = x3_2 + math.tan(alpha)*l3

p1 = {x1_1,y1}
p2 = {x2_1,y2}
p3 = {x2_2 ,y2}
p4 = {x3_3,y3+1} 
p5 = {x3_2,y3}
p6 = {x4_2,y4}
p7 = {x4_1,y4}
p8 = {x3_1,y3}

-- this function produces negative value
 
function mx(p)
   return { - p[1], p[2] }
end

points1 = { p1, p2, p3, p4, mx(p4), mx(p3), mx(p2), mx(p1) }
points2 = { p5,p6,p7,p8 }
points3 = { mx(p5),mx(p6),mx(p7),mx(p8) }

--points = { p1, p2, p3, p4, p5, p6, p7, p8, mx(p8), mx(p7), mx(p6), mx(p5), mx(p4), mx(p3), mx(p2), mx(p1) }

for i,v in ipairs(points1) do
   print("p["..tostring(i).."] = ",v[1],v[2])
end

joffs = height_bsio2+height_bsi+height_wg/2.

prism1 = Transform{ 
   Prism{  
      points = points1,
      height = height_wg
   }, 
   axis = { 1, 0, 0 },
   angle = 90,
   move = { 0,joffs,0 } 
}

prism2 = Transform{ 
   Prism{  
      points = points2,
      height = height_wg
   },
   axis = { 1, 0, 0 },
   angle = 90,
   move = { 0,joffs,0 } 
}
prism3 = Transform{ 
   Prism{  
      points = points3,
      height = height_wg
   },
   axis = { 1, 0, 0 },
   angle = 90,
   move = { 0,joffs,0 } 
}

mmi:add{ prism1, depth = 1, value=eps_si }
mmi:add{ prism2, depth = 1, value=eps_si }
mmi:add{ prism3, depth = 1, value=eps_si }
