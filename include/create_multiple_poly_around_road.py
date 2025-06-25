import numpy as np

def create_polys_around_road(xa,ya,xb,yb,width):
    from include.create_poly_around_road import create_2_poly_around_road
    from shapely.geometry import LineString
    road = LineString(((xa, ya), (xb, yb)))
    # To get the line length in meters, project it in Pseudo-Mercator
    import shapely.ops as sp_ops
    import pyproj
    my_transformer = pyproj.Transformer.from_crs('EPSG:2154', 'EPSG:3857', always_xy=True) # from lamb93 to Pseudo-Mercator, NB : latlon='EPSG:4326'
    line_transformed = sp_ops.transform(my_transformer.transform, road)
    line_length = line_transformed.length
    #line_length = np.sqrt((xb - xa)**2 + (yb-ya)**2)
    # Depending on the line length, define a given number of subpolygons
    road_polys = []
    if line_length < 25 :
        n_poly = 1  # one polygon for each side of the road
        road_polys.extend(create_2_poly_around_road(xa,ya,xb,yb,width))
    ## if the road is longer than 10m, create several polygons
    else :
        if 25 <= line_length < 50:
            n_poly = 5  # 5 subparts of the road, with two sides
        else :
            n_poly = 10  # 10 subparts of the road, with two sides
        dist = float(line_length / n_poly)
        for i in range (1, n_poly+1):  # for (n-1) subpoly, set C between A and B and create the subpoly
            if i < n_poly :
                (xc, yc) = road.interpolate(dist*i).coords[0]
            else :  # for the nth subpoly, create the subpoly between the (n-1)th C and B
                (xc, yc) = (xb, yb)
            road_polys.extend(create_2_poly_around_road(xa, ya, xc, yc, width))
            (xa,ya) = (xc,yc)
    n_poly = 2*n_poly
    return([road_polys,n_poly,line_length])

def create_left_right_polys_around_road(xa,ya,xb,yb,width):
    from include.create_poly_around_road import create_2_poly_around_road
    from shapely.geometry import LineString, Point
    from shapely import centroid
    import shapely.ops as sp_ops
    import pyproj

    road = LineString(((xa, ya), (xb, yb)))
    # To get the line length in meters, project it in Pseudo-Mercator
    my_transformer = pyproj.Transformer.from_crs('EPSG:2154', 'EPSG:3857', always_xy=True) # from lamb93 to Pseudo-Mercator, NB : latlon='EPSG:4326'
    line_transformed = sp_ops.transform(my_transformer.transform, road)
    line_length = line_transformed.length
    #line_length = np.sqrt((xb - xa)**2 + (yb-ya)**2)
    # Depending on the line length, define a given number of subpolygons
    road_polys_left = []
    road_polys_right = []
    centers = []
    if line_length < 15 :
        n_poly = 1  # one polygon for each side of the road
        [left_poly, right_poly] = create_2_poly_around_road(xa,ya,xb,yb,width)
        road_polys_left.append(left_poly) # extend
        road_polys_right.append(right_poly) # extend
        centers.append(centroid(road))
    ## if the road is longer than 10m, create several polygons
    else :
        if line_length < 150:
            n_poly = int(np.rint(line_length/10)) # 10m long subparts of the road, with two sides
        else :
            n_poly = 10  # 10 subparts of the road, with two sides
        dist = float(line_length / n_poly)
        for i in range (1, n_poly+1):  # for (n-1) subpoly, set C between A and B and create the subpoly
            if i < n_poly :
                (xc, yc) = road.interpolate(dist*i).coords[0]
            else :  # for the nth subpoly, create the subpoly between the (n-1)th C and B
                (xc, yc) = (xb, yb)

            subroad = LineString(((xa, ya), (xc, yc)))
            centers.append(centroid(subroad))# centers.append(Point(xc, yc))
            [left_poly, right_poly] = create_2_poly_around_road(xa, ya, xc, yc, width)
            road_polys_left.append(left_poly)  # extend
            road_polys_right.append(right_poly)  # extend
            (xa,ya) = (xc,yc)
    return([road_polys_left, road_polys_right,centers, n_poly, line_length])
