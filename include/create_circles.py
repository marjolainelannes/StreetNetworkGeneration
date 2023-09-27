def create_circles_around_road_nodes(xa, xb, ya, yb, radius) :
    from shapely.geometry import Point, Polygon
    node_A = Point(xa,ya)
    node_B = Point(xb, yb)
    circle_A = node_A.buffer(radius) # type is polygon
    circle_B = node_B.buffer(radius)
    return(circle_A, circle_B)

