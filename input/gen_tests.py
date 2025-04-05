import gen_circle_regions as gcirc
import gen_rectangle_regions as grect
import os

SCENE = 10
VERTICES = 9
MIN_SIDE = 20
MAX_SIDE = 40
MIN_RADIOUS = 10
MAX_RADIOUS = 20
N_MESH_POINTS = 1000
N_INSTANCES = 10
DIR_PATH = "tests/"

os.makedirs(DIR_PATH, exist_ok=True)

for i in range(N_INSTANCES):
    # Circle instance
    circles = gcirc.gen_circles(
        SCENE, VERTICES, MIN_RADIOUS, MAX_RADIOUS, real=True, seed=i
    )
    points = gcirc.gen_mesh_points(SCENE, N_MESH_POINTS)  # 1M points
    num_vertices, num_edges, edges = gcirc.build_hypergraph(circles, points)
    name = f"circ_{SCENE}_{VERTICES}_{MIN_RADIOUS}_{MAX_RADIOUS}_{i}"
    gcirc.write_hypergraph(
        DIR_PATH + name + ".txt",
        (num_vertices, num_edges, edges),
    )
    gcirc.plot(SCENE, circles, saveto=DIR_PATH + name + ".png")
    # Rectangular instance
    rectangles = grect.gen_rectangles(
        SCENE, VERTICES, MIN_SIDE, MAX_SIDE, real=True, seed=i
    )
    num_vertices, num_edges, edges = grect.build_hypergraph(rectangles)
    grect.write_hypergraph(
        DIR_PATH + f"rect_{SCENE}_{VERTICES}_{MIN_SIDE}_{MAX_SIDE}_{i}.txt",
        (num_vertices, num_edges, edges),
    )
