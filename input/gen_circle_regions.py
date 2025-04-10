"""Random generation of hypergraphs associated with circle regions."""

import math
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Type alias
Number = float | int
Point = tuple[Number, Number]
Circ = tuple[Point, Number]  # x-coordinate, y-coordinate, radious
Edge = list[int]
HyperGraph = tuple[int, int, tuple[int, list[int]]]
Vertex = int


def gen_centers(
    side: int,
    ncenters: int,
    real: bool = False,
    seed: int | None = None,
) -> list[Point]:
    """
    Generate a list of centers.

    Args:
        side: Length of the side of the scene.
        ncenters: Number of circles.
        real: Whether to generate random real numbers, otherwise integers.
        seed: Seed for random generation.

    Returns:
        A list of centers. Each one of the form (x,y) where (x,y) are
        the coordinates of the center.
    """
    ncirc = 0
    circ = []
    generator = random.uniform if real else random.randint
    random.seed(seed)
    while ncirc < ncenters:
        x = generator(0, side)
        y = generator(0, side)
        circ.append((x, y))
        ncirc += 1
    return circ


def gen_circles(
    side: int,
    centers: list[Point],
    radius: int,
) -> list[Circ]:
    """
    Generate a list of circles.

    Args:
        side: Length of the side of the scene.
        centers: List of centers.
        radius: Length of the radius,
            written as a percentage of side.

    Returns:
        A list of circles. Each one of the form (x,y,r) where (x,y) are
        the coordinates of the center and r is the radious.
    """
    return [(c, radius * side / 100) for c in centers]


def plot(side: int, circ: list[Circ], saveto: str | None = str) -> None:
    """Plot a list of circles inside a scene.

    Args:
        side: Length of the side of the scene.
        circ: List of rectangles.
        saveto: Path to save plot.
    """
    plt.figure()
    plt.xlim(0, side)
    plt.ylim(0, side)
    ax = plt.gca()
    for c, r in circ:
        ax.add_patch(Circle(c, r, fill=False))
    if saveto is None:
        plt.show()
    else:
        plt.savefig(saveto)


def gen_mesh_points(side: int, nsteps: int) -> list[Point]:
    """Generate the points of a mesh.

    Args:
        side: Length of the side of the scene.
        nsteps: Size of the mesh, that is,
            with (nsteps + 1)*(nsteps + 1) points.
    """
    step = side / nsteps
    return [(x * step, y * step) for x in range(nsteps + 1) for y in range(nsteps + 1)]


def build_hypergraph(circ: list[Circ], points: list[Point]) -> HyperGraph:
    """Build the hypergraph associated with the circle regions.

    Args:
        circ: A list of circles.
        points: A list of points in R^2.

    Returns:
        A hypergraph. That is, a tuple with: number of vertices, number of
        hyperedges, and a list of hyperedges. Each hyperedge is a tuple with
        the number of implied vertices and a list of them.
    """
    n = len(circ)
    m = 0
    edges = []
    for p in points:
        nedge = 0
        edge = []
        for i, (c, r) in enumerate(circ):
            if euclidean_distance(p, c) > r:
                continue
            nedge += 1
            edge.append(i)
        if nedge > 0 and (nedge, edge) not in edges:
            m += 1
            edges.append((nedge, edge))
    return n, m, edges


def write_hypergraph(path: str, graph: HyperGraph) -> None:
    """Write a hypergraph in a file.

    Args:
        path: Path of file.
        graph: Hypergraph.
    """
    with open(path, "w", encoding="utf8") as f:
        f.write(str(graph[0]) + " " + str(graph[1]) + "\n")
        for e in graph[2]:
            f.write(str(e[0]))
            for v in e[1]:
                f.write(" " + str(v))
            f.write("\n")


def euclidean_distance(point1: Point, point2: Point) -> float:
    """Euclidean distance between two points in R2.

    Args:
        point1: Point in R2.
        point2: Point in R2.

    Returns:
        Euclidean distance.
    """
    return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
