"""Random generation of hypergraphs associated with axis-parallel rectangle
regions."""

import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Type alias
Number = float | int
Rect = tuple[Number, Number, Number, Number]
Edge = list[int]
HyperGraph = tuple[int, int, tuple[int, list[int]]]


def gen_rectangles(
    sside: int,
    rnum: int,
    rmin: int,
    rmax: int,
    real: bool = False,
    seed: int | None = None,
) -> list[Rect]:
    """
    Generate a list of random rectangles.

    Args:
        sside: Length of the side of the scene.
        rnum: Number of rectangles.
        rmin: Minimum length of the side of the rectangles,
            written as a percentage of sside.
        rmax: Maximum length of the side of the rectangles
            written as a percentage of sside.
        real: Whether to generate random real numbers, otherwise integers.
        seed: Seed for random generation.

    Returns:
        A list of rectangles. Each one of the form (x1,y1,x2,y2) where (x1,y1)
        and (x2,y2) are the coordinates of the opposite vertices.
    """
    nrect = 0
    rect = []
    generator = random.uniform if real else random.randint
    random.seed(seed)
    while nrect < rnum:
        x1 = generator(0, sside)
        y1 = generator(0, sside)
        x2 = x1 + generator(sside * rmin / 100, sside * rmax / 100)
        y2 = y1 + generator(sside * rmin / 100, sside * rmax / 100)
        if x2 > sside or y2 > sside:
            continue
        rect.append((x1, y1, x2, y2))
        nrect += 1
    return rect


def plot(sside: int, rect: list[Rect]) -> None:
    """Plot a list of rectangles inside a scene.

    Args:
        sside: Length of the side of the scene.
        rect: List of rectangles.
    """
    plt.figure()
    plt.xlim(0, sside)
    plt.ylim(0, sside)
    ax = plt.gca()
    for x1, y1, x2, y2 in rect:
        ax.add_patch(Rectangle((x1, y1), x2 - x1, y2 - y1, fill=False))
    plt.show()


def get_axis_coordinates(rect: list[Rect], axis: str) -> list[Number]:
    """
    Given a axis, get all coordinates where a rectangle begins or ends.

    Args:
        rect: List of rectangles.
        axis: Name of axis: "x" or "y".

    Returns:
        A list of coordinates.
    """
    if axis == "x":
        return [x for r in rect for x in (r[0], r[2])]
    return [y for r in rect for y in (r[1], r[3])]


def get_midpoints(numbers: list[Number]) -> list[float]:
    """Get the midpoints between every pair of consecutive elements of a list.

    Args:
        l: A sorted list of numbers.

    Returns:
        A list of midpoints.
    """
    return list(map(lambda x, y: (x + y) / 2, numbers[:-1], numbers[1:]))


def build_hypergraph(rect: list[Rect]) -> HyperGraph:
    """Generate the hypergraph associated with the rectangle regions.

    Args:
        rect: A list of rectangles.

    Returns:
        A hypergraph. That is, a tuple with: number of vertices, number of
        hyperedges, and a list of hyperedges. Each hyperedge is a tuple with the
        number of implied vertices and a list of them.
    """
    n = len(rect)
    m = 0
    edges = []
    mx = get_midpoints(list(set(get_axis_coordinates(rect, "x"))))
    my = get_midpoints(list(set(get_axis_coordinates(rect, "y"))))
    for x in mx:
        for y in my:
            nedge = 0
            edge = []
            for i, (x1, y1, x2, y2) in enumerate(rect):
                if x < x1 or x > x2 or y < y1 or y > y2:
                    continue
                nedge += 1
                edge.append(i)
            if nedge > 0 and (nedge, edge) not in edges:
                m += 1
                edges.append((nedge, edge))
    return n, m, edges


def write_hypergraph(path: str, graph: HyperGraph):
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
