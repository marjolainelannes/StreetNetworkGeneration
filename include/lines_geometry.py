def create_network_graph(links_df, nodes_df, crs="", keep_directed = True):
    import networkx as nx
    # Create graph
    G_network = nx.from_pandas_edgelist(
        links_df, source="from", target="to", edge_attr=True, create_using=nx.MultiDiGraph()
    )
    if len(crs) > 0 :
        G_network.graph["crs"] = crs
    nodes_df["id"] = nodes_df["id"].astype(str)
    G_network.add_nodes_from(nodes_df["id"])
    # Add node attributes to the graph
    for node in nodes_df.itertuples():
        G_network.nodes[node.id]["x"] = node.x
        G_network.nodes[node.id]["y"] = node.y
    # undirected option
    if not keep_directed :
        G_network = G_network.to_undirected()
    return(G_network)

def project_point_on_line(point_to_project : list, line_pointA : list, line_pointB : list) :
    from skspatial.objects import Line
    line = Line.from_points(line_pointA, line_pointB)
    projected_point_coordinates = line.project_point(point_to_project)
    return(projected_point_coordinates)

def distance_to_line(point_coordinates : list, line_pointA : list, line_pointB : list) :
    from skspatial.objects import Line
    line_1 = Line.from_points(line_pointA, line_pointB)
    distance = line_1.distance_point(point_coordinates)
    #print(line_1.project_point(point_coordinates))
    return(distance)

def lines_intersection(line1_pointA : list, line1_pointB : list, line2_pointC : list, line2_pointD : list):
    # returns an object Point from scikit spatial
    from skspatial.objects import Line
    line_1 = Line.from_points(line1_pointA, line1_pointB)
    line_2 = Line.from_points(line2_pointC, line2_pointD)
    intersection_coordinates = line_2.intersect_line(line_1)
    return(intersection_coordinates)

def slope(coordinates:list): # Line slope given two points:
    [x1, y1, x2, y2] = coordinates
    s = (y2-y1)/(x2-x1)
    return(s)

def angle_two_lines(slope1, slope2):
    import math
    angle = math.degrees(math.atan((slope2-slope1)/(1+(slope2*slope1))))
    return(angle)

def segments_on_the_same_line(line1_pointA : list, line1_pointB : list, line2_pointC : list, line2_pointD : list) :
    import math
    line1 = line1_pointA + line1_pointB
    line2 = line2_pointC + line2_pointD
    s1 = slope(line1)
    s2 = slope(line2)
    theta = angle_two_lines(s1, s2)
    dist = distance_to_line(line1_pointA, line2_pointC, line2_pointD)
    #print("distance", dist, "angle", theta, "sinus", abs(math.sin(math.radians(theta))))
    if abs(math.sin(math.radians(theta))) < 0.5 : # NB: sin(30) = 0.5
        if (abs(math.sin(math.radians(theta))) > 0.1736)  and (dist < 10): # NB: sin(10) = 0.1736
            return(True)
        # flat angles require a smaller distance to check if the lines are not parallel
        elif (abs(math.sin(math.radians(theta))) < 0.1736) and (dist < 5) :
            return(True)
        else :
            return(False)
    else :
        return(False)

def parallel_lines(points_coordinates_segment_1 : list, points_coordinates_segment_2 : list) :
    import math
    s1 = slope(points_coordinates_segment_1)
    s2 = slope(points_coordinates_segment_2)
    theta = angle_two_lines(s1, s2)
    if abs(math.sin(math.radians(theta))) < 0.5 : # NB: sin(30) = 0.5
        return(True)
    else :
        return(False)

def get_street_line(node : list, lines_A_list : list, lines_B_list : list) :
    # node : [x,y] to project on all lines
    # lines_A_list : list of coordinates [xa, ya] of one point A of each line
    # lines_B_list : list of coordinates [xb, yb] of one point B of each line
    import numpy as np
    from skspatial.objects import Points
    # 1. Calculate the slope of the street : mean slope of lines
    n_lines = len(lines_A_list)
    slopes = []
    for k in range(n_lines):
        slopes.append(slope([lines_A_list[k][0], lines_A_list[k][1], lines_B_list[k][0], lines_B_list[k][1]]))
    a = np.mean(slopes)
    # 2. Project the node on all the lines and return the centroid
    projected_points_list = []
    for i in range(n_lines) :
        new_point = project_point_on_line(node, lines_A_list[i], lines_B_list[i])
        projected_points_list.append(new_point)
    projected_points = Points(projected_points_list)
    points_centered, centroid = projected_points.mean_center(return_centroid=True)
    # 3. calcule the y_intercept of the line
    y_intercept = centroid[1] - a * centroid[0]
    return(a, y_intercept)

# a tester : get_street_line
