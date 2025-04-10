import gen_circle_regions as gcirc
import os

SIDE = 1
VERTICES = 10
RADIOUS = [40]
N_MESH_POINTS = 1000
N_INSTANCES = 10
DIR_PATH = "tests/"

os.makedirs(DIR_PATH, exist_ok=True)

for i in range(N_INSTANCES):

    # Circle instances
    centers = gcirc.gen_centers(SIDE, VERTICES, real=True, seed=i)
    points = gcirc.gen_mesh_points(SIDE, N_MESH_POINTS)  # 1M points
    for r in RADIOUS:
        circles = gcirc.gen_circles(SIDE, centers, r)
        num_vertices, num_edges, edges = gcirc.build_hypergraph(circles, points)
        NAME = f"circ_{SIDE}_{VERTICES}_{r}_{i}"
        gcirc.write_hypergraph(
            DIR_PATH + NAME + ".txt",
            (num_vertices, num_edges, edges),
        )
        gcirc.plot(SIDE, circles, saveto=DIR_PATH + NAME + ".png")
