module ("funcs", package.seeall)


function roundbox(hlength,hwidth,hheight,redges)

  local rndbx = {}
  
  box1 = Box{
    from = {-hlength,-hwidth,-hheight+redges},
    to = {hlength,hwidth,hheight-redges}
  }
  box2 = Box{
    from = {-hlength,-hwidth+redges,-hheight},
    to = {hlength,hwidth-redges,hheight}
  }
  box3 = Box{
    from = {-hlength+redges,-hwidth,-hheight},
    to = {hlength-redges,hwidth,hheight}
  }
  rndbx = BinaryOr{box1,box2}
  rndbx = BinaryOr{rndbx,box3}

  -- round edges by removing box edges and adding cylinders
  for i = -1,1,2 do
    for j = -1,1,2 do
      for k = -1,1,2 do
        sphere = Sphere{
          at = {i*(hlength-redges),j*(hwidth-redges),k*(hheight-redges)},
          radius = {redges}
        }
        rndbx = BinaryOr{rndbx,sphere}
      end
    end
  end

  for i = -1,1,2 do
    for j = -1,1,2 do
      cyl1 = Cylinder{
        at = {0,i*(hwidth-redges),j*(hheight-redges)},
        radius = redges,
        hheight = 2*(hlength-redges)
      }
      rndbx = BinaryOr{rndbx,cyl1}
      cyl1 = Cylinder{
        at = {i*(hlength-redges),0,j*(hheight-redges)},
        radius = redges,
        hheight = 2*(hwidth-redges)
      }
      rndbx = BinaryOr{rndbx,cyl1}
      cyl1 = Cylinder{
        at = {i*(hlength-redges),j*(hwidth-redges),0},
        radius = redges,
        hheight = 2*(hheight-redges)
      }
      rndbx = BinaryOr{rndbx,cyl1}

    end
  end
  
  return rndbx

end

