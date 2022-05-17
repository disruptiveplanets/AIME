import numpy as np

NUM_DRAWS = 10
SHIFT_FRACTION = 1

def get_grids_centroid(num_grids, grid_line, indicator_map, indicator):
    # Shift each point to its cell's centroid and iterate
    points = get_rand_points(num_grids, grid_line, indicator)
    x,y,z = np.meshgrid(grid_line, grid_line, grid_line)
    grids = None
    for i in range(NUM_DRAWS):
        grids = get_grids_voronoi(points, grid_line, indicator_map)
        for grid_index, grid in enumerate(grids):
            centroid = np.array([np.sum(grid * x), np.sum(grid * y), np.sum(grid * z)]) / np.sum(grid)
            points[grid_index] = points[grid_index] + (centroid - points[grid_index]) * SHIFT_FRACTION
        #print(np.sum(grids, axis=(1,2,3)))
    grids = get_grids_voronoi(points, grid_line, indicator_map)
    return points, grids

def get_grids_weights(num_grids, grid_line, indicator, indicator_map):
    # Weight the points with small volumes stronger and iterate
    points = get_rand_points(num_grids, grid_line, indicator)
    weights = np.ones(num_grids)
    grids = None
    equal_vol = np.sum(indicator_map) / num_grids
    for i in range(NUM_DRAWS):
        grids = get_grids_weight(points, weights)
        vols = np.sum(grids, axis=(1,2,3)) / equal_vol
        weights = (1 / vols)**0.3
        print(vols)
    grids = get_grids_weight(points, weights)
    return points, grids

def get_rand_points(num_grids, grid_line, indicator):
    points = []
    while len(points) < num_grids:
        p = np.random.uniform(low=np.min(grid_line), high=np.max(grid_line), size=3)
        if indicator(p[0], p[1], p[2]):
            points.append(p)
    return np.array(points)

def get_grids_voronoi(points, grid_line, indicator_map):
    x,y,z = np.meshgrid(grid_line, grid_line, grid_line)
    dists = np.zeros((len(points), len(grid_line), len(grid_line), len(grid_line)))
    for i in range(len(points)):
        dists[i] = (x - points[i,0]) ** 2 + (y - points[i,1]) ** 2 + (z - points[i,2]) ** 2
    grid_indices = np.argmin(dists, axis=0)

    grids = []
    for i in range(len(points)):
        grids.append((grid_indices == i) * indicator_map)

    return np.array(grids)

def get_grids_weight(points, weights, grid_line, indicator_map):
    x,y,z = np.meshgrid(grid_line, grid_line, grid_line)
    dists = np.zeros((len(points), len(grid_line), len(grid_line), len(grid_line)))
    for i in range(len(points)):
        dists[i] = (x - points[i,0]) ** 2 + (y - points[i,1]) ** 2 + (z - points[i,2]) ** 2 / weights[i]**2
    grid_indices = np.argmin(dists, axis=0)

    grids = []
    for i in range(len(points)):
        grids.append((grid_indices == i) * indicator_map)

    return np.array(grids)