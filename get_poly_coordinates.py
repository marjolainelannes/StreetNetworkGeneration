# Function to get the coordinates of the 4 summits of the polygon surrounding the road with xa =/ xb and ya =/ yb
# Function to get a point on a line at a given distance from A
def xycoord_dist_x0_y0(a,xa,ya,xb,yb,dist):
    #x,y coordinates on line2 with distance dist from (x0,y0)
    import numpy as np
    if a > 0 :
       x1 = xa - dist * np.sqrt(a**2/(a**2+1))
       y1 = ya + dist / (np.sqrt(a**2+1))

       x2 = xa + dist * np.sqrt(a**2/(a**2+1))
       y2 = ya - dist / (np.sqrt(a**2+1))

       x3 = xb + dist * np.sqrt(a**2/(a**2+1))
       y3 = yb - dist / (np.sqrt(a**2+1))

       x4 = xb - dist * np.sqrt(a**2/(a**2+1))
       y4 = yb + dist / (np.sqrt(a**2+1))
    if a < 0 :
       x1 = xa - dist * np.sqrt(a**2/(a**2+1))
       y1 = ya - dist / (np.sqrt(a**2+1))

       x2 = xa + dist * np.sqrt(a**2/(a**2+1))
       y2 = ya + dist / (np.sqrt(a**2+1))

       x3 = xb + dist * np.sqrt(a**2/(a**2+1))
       y3 = yb + dist / (np.sqrt(a**2+1))

       x4 = xb - dist * np.sqrt(a**2/(a**2+1))
       y4 = yb - dist / (np.sqrt(a**2+1))
    return ((x1,y1),(x2,y2),(x3,y3),(x4,y4))

# à supprimer : il y a déjà la fonction interpolate qui le fair
def get_C_from_line_dist(line, dist) :
    from shapely.geometry import LineString  # usefull for line.coords
    (xc, yc) = line.interpolate(dist).coords[0]
    # slope = (yb - ya) / (xb - xa)
    # import numpy as np
    # num_vert = int(round(line.length / dist))
    # if num_vert == 0:
    #     num_vert = 1
    #     # geom.interpolate(dist)
    # return LineString(
    #     [line.interpolate(float(n) / num_vert, normalized=True)
    #      for n in range(num_vert + 1)])
    return([xc,yc])