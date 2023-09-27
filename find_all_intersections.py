def find_all_intersections(roads_bdtopo_df, poly):
    from shapely.geometry import Polygon, Point, LineString
    ids_list = []
    intersections_list = []
    widths_list = []
    lengths_list = []
    for i, road_geom in enumerate(roads_bdtopo_df['LINE']):
        if poly.intersects(road_geom):
            id = roads_bdtopo_df.loc[i, 'ROAD_ID']
            width = roads_bdtopo_df.loc[i, 'WIDTH']
            road_length = road_geom.length
            lengths_list.append(road_length)
            ids_list.append(id)
            intersections_list.append(poly.intersection(road_geom).length)
            widths_list.append(width)
    return(ids_list, intersections_list, lengths_list,widths_list)