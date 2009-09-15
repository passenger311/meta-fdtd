module ("funcs", package.seeall)


function roundbox(hlength,hwidth,hheight,redges)

  local rndbx = {}
  
  box1 = Box{
    from = {-hlength,-hwidth+redges,-hheight+redges},
    to = {hlength,hwidth-redges,hheight-redges}
  }
  box2 = Box{
    from = {-hlength+redges,-hwidth+redges,-hheight},
    to = {hlength-redges,hwidth-redges,hheight}
  }
  box3 = Box{
    from = {-hlength+redges,-hwidth,-hheight+redges},
    to = {hlength-redges,hwidth,hheight-redges}
  }
  rndbx = BinaryOr{box1,box2}
  rndbx = BinaryOr{rndbx,box3}

  -- round edges by removing box edges and adding cylinders
  for i = -1,1,2 do
    for j = -1,1,2 do
      for k = -1,1,2 do
        sphere = Sphere{
          at = {i*(hlength-redges),j*(hwidth-redges),k*(hheight-redges)},
          radius = redges
        }
        rndbx = BinaryOr{rndbx,sphere}
      end
    end
  end

  for i = -1,1,2 do
    for j = -1,1,2 do
      cyl1 = Cylinder{
        at = {0,0,0},
        radius = redges,
        height = 2*(hlength-redges)
      }
      cyl1 = Transform{
        cyl1,
        move = {0,i*(hwidth-redges),j*(hheight-redges)},
        angle = 90,
        axis = {0,1,0}
      }
      rndbx = BinaryOr{rndbx,cyl1}

      cyl2 = Cylinder{
        at = {0,0,0},
        radius = redges,
        height = 2*(hwidth-redges)
      }
      cyl2 = Transform{
        cyl2,
        move = {i*(hlength-redges),0,j*(hheight-redges)},
        angle = 90,
        axis = {1,0,0}
      }
      rndbx = BinaryOr{rndbx,cyl2}

      cyl3 = Cylinder{
        at = {0,0,0},
        radius = redges,
        height = 2*(hheight-redges)
      }
      cyl3 = Transform{
        cyl3,
        move = {i*(hlength-redges),j*(hwidth-redges),0},
      }
      rndbx = BinaryOr{rndbx,cyl3}

    end
  end

  return rndbx

end

