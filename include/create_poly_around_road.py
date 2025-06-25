from shapely.geometry import Polygon

# Function to get the coordinates of the 4 summits of the polygon surrounding the road with xa =/ xb and ya =/ yb
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

def create_2_poly_around_road(xa,ya,xb,yb,width):
    # Get the polygon coordinates
    if xa == xb :  # vertical axis
                (x1,y1)=(xa-width,ya)
                (x2,y2)=(xa+width,ya)
                (x3,y3)=(xa+width,yb)
                (x4,y4)=(xa-width,yb)
    elif ya == yb :  # horizontal axis
                (x1,y1)=(xa,ya-width)
                (x2,y2)=(xa,ya+width)
                (x3,y3)=(xb,ya+width)
                (x4,y4)=(xb,ya-width)
    else:
                slope=(yb-ya)/(xb-xa)
                (x1, y1), (x2, y2), (x3, y3), (x4, y4) = xycoord_dist_x0_y0(slope, xa, ya, xb, yb, width)
    # Create one polygon for each side of the road
    road_poly_1 = Polygon(((x1, y1), (xa, ya), (xb, yb), (x4, y4)))
    road_poly_2 = Polygon(((xa, ya), (x2, y2), (x3, y3), (xb, yb)))
    return([road_poly_1, road_poly_2])

def create_1_poly_around_road(xa,ya,xb,yb,width):
    # Get the polygon coordinates
    if xa == xb :  # vertical axis
                (x1,y1)=(xa-width,ya)
                (x2,y2)=(xa+width,ya)
                (x3,y3)=(xa+width,yb)
                (x4,y4)=(xa-width,yb)
    elif ya == yb :  # horizontal axis
                (x1,y1)=(xa,ya-width)
                (x2,y2)=(xa,ya+width)
                (x3,y3)=(xb,ya+width)
                (x4,y4)=(xb,ya-width)
    else:
                slope=(yb-ya)/(xb-xa)
                (x1, y1), (x2, y2), (x3, y3), (x4, y4) = xycoord_dist_x0_y0(slope, xa, ya, xb, yb, width)
    # Create one polygon
    road_poly=Polygon(((x1, y1), (x2, y2), (x3, y3), (x4, y4)))
    return(road_poly)
    
