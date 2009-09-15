module ("funcs", package.seeall)


function roundstrip(hwidth,hheight,redges)

  local rndstrp = {}

  strip = Box{
    from = {-hwidth,-hheight,-2},
    to = {hwidth,hheight,2}
  }

  -- round edges by removing box edges and adding cylinders

  strip1 = Box{
    from = {-hwidth,-hheight,-2},
    to = {-hwidth+redges,-hheight+redges,2}
  }
  strip = BinaryAndNot{strip,strip1}
  strip1 = Box{
    from = {hwidth-redges,-hheight,-2},
    to = {hwidth,-hheight+redges,2}
  }
  strip = BinaryAndNot{strip,strip1}
  strip1 = Box{
    from = {-hwidth,hheight-redges,-2},
    to = {-hwidth+redges,hheight,2}
  }
  strip = BinaryAndNot{strip,strip1}
  strip1 = Box{
    from = {hwidth-redges,hheight-redges,-2},
    to = {hwidth,hheight,2}
  }
  strip = BinaryAndNot{strip,strip1}

  cyl = Cylinder{
    at = {-hwidth+redges,-hheight+redges,0},
    radius = redges,
    height = 4
  }
  strip = BinaryOr{strip,cyl}
  cyl = Cylinder{
    at = {-hwidth+redges,hheight-redges,0},
    radius = redges,
    height = 4
  }
  strip = BinaryOr{strip,cyl}
  cyl = Cylinder{
    at = {hwidth-redges,-hheight+redges,0},
    radius = redges,
    height = 4
  }
  strip = BinaryOr{strip,cyl}
  cyl = Cylinder{
    at = {hwidth-redges,hheight-redges,0},
    radius = redges,
    height = 4
  }
  rndstrp = BinaryOr{strip,cyl}
  -- finished rounding edges
  return rndstrp

end