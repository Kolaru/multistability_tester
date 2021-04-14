import networkx as nx
import numpy as np

from itertools import combinations
from matplotlib import pyplot as plt
from matplotlib.patches import ArrowStyle

from utils import pairs

plt.ioff()

class MutableEdge:
    def __init__(self, ax,
                 nedges=0, k=0,
                 name="", sign=0,
                 v1=None, v2=None,
                 patchA=None, patchB=None,
                 sep=0):

        self.v1 = v1
        self.v2 = v2
        self.k = k
        self.sign = sign
        self.patchA = patchA
        self.patchB = patchB
        self.ax = ax

        self._sep = sep

        self.dv = v2 - v1
        self.cross_dv = np.array([self.dv[1], -self.dv[0]])  # Perpendicular to self.dv
        self.mid = 0.5*(v1 + v2)
        self.norm_dv = np.sqrt(self.dv[0]**2 + self.dv[1]**2)
        self.arrows = []

        self.add_arrow()
        self.label = ax.text(*self.label_pos, name, va="center", ha="center",
                             bbox=dict(facecolor="w",
                                       edgecolor="none"),
                             zorder=100)
    @property
    def rad(self):
        return 2*self.k*self.sep/self.norm_dv

    @property
    def label_pos(self):
        return self.mid + 0.5*self.cross_dv*self.rad

    @property
    def connection_style(self):
        return "arc3, rad={}".format(self.rad)

    @property
    def sep(self):
        return self._sep

    @sep.setter
    def sep(self, sep):
        self._sep = sep
        for arrow in self.arrows:
            arrow.set_connectionstyle(self.connection_style)
        self.label.set_position(self.label_pos)

    def add_arrow(self, highlight=False):
        if highlight:
            color = "yellow"
            arrowstyle = "-"
            linewidth = 8
        else:
            linewidth = 2
            if self.sign > 0:
                color = "green"
                arrowstyle = "->"
            elif self.sign < 0:
                color = "red"
                arrowstyle= "|-|, widthA=0, widthB=0.5"
            else:
                color = "orange"
                arrowstyle = "-"

        arrowprops = dict(arrowstyle=arrowstyle,
                          edgecolor=color,
                          linewidth=linewidth,
                          shrinkA=0, shrinkB=5,
                          patchA=self.patchA,
                          patchB=self.patchB,
                          connectionstyle=self.connection_style)

        edge = self.ax.annotate("",
                                xy=self.v2,
                                xytext=self.v1,
                                arrowprops=arrowprops)

        if highlight:
            edge.set_zorder(-10000)

        self.arrows.append(edge.arrow_patch)


class MutableInfluenceGraph:
    def __init__(self, ax, edges, separation):
        self.figure = ax.figure
        self.ax = ax
        self.edges = edges
        self.separation = separation

    def highlight_hooping(self, cycles, reactions):
        for cycle, path in zip(cycles, reactions):
            for e, R in zip(pairs(cycle), path):
                mutable = self.edges[e][R]
                mutable.add_arrow(highlight=True)

    def update_sep(self, sep):
        for e in self.edges.values():
            for arrow in e.values():
                arrow.sep = sep*self.separation

        self.update()

    def update(self):
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

def plot_influence_graph(influence_graph):
    # Create the figure
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.axis('off')  # Remove axes
    ax.set_aspect('equal')  # Avoid deforming the graph

    layout = nx.kamada_kawai_layout(influence_graph)  # This determines where the node will go

    # Find nodes that are the furthest away and rotate to make this axis horizontal
    vectors = ((v1, v2) for v1, v2 in combinations(layout.values(), 2))
    v1, v2 = max(vectors, key=lambda vecs: np.sum((vecs[1]-vecs[0])**2))
    v = v2 - v1

    edges_length = (np.sum((layout[a] - layout[b])**2) for (a, b) in influence_graph.edges())
    max_edge_length = max(edges_length)
    sep = 0.1 * max_edge_length

    angle = np.arctan2(v[1], v[0])
    rot = np.array([[ np.cos(angle), np.sin(angle)],
                    [-np.sin(angle), np.cos(angle)]])

    layout = {key:rot @ val for key, val in layout.items()}

    xs = [v[0] for v in layout.values()]
    ys = [v[1] for v in layout.values()]

    minx = min(xs)
    maxx = max(xs)

    miny = min(ys)
    maxy = max(ys)

    pad = 0.05*max(maxx - minx, maxy - miny)

    ax.set_xlim(minx - pad, maxx + pad)
    ax.set_ylim(miny - pad, maxy + pad)

    node_patches = {}
    biggest_label = max(influence_graph.nodes, key=lambda x:len(x))
    txtstyle = dict(ha="center", va="center", fontsize=12)

    for label, pos in layout.items():
        txt = ax.text(pos[0], pos[1], biggest_label, color="none",
                      bbox=dict(boxstyle="circle, pad=0.5", alpha=1, facecolor="white"),
                      **txtstyle)

        patch = txt.get_bbox_patch()
        node_patches[label] = patch

        ax.text(pos[0], pos[1], label, **txtstyle)

    nx.draw_networkx_nodes(influence_graph, layout, ax=ax, node_color="none")

    plotted_pairs = []
    edges = {}

    for a, b, data in influence_graph.edges(data=True):
        if (b, a) in plotted_pairs:
            continue

        edges[(a, b)] = {}
        edges[(b, a)] = {}
        reactions = data["reactions"]
        sign = data["sign"]

        if a == b:
            if sign > 0:
                mec = "green"
            elif sign < 0:
                mec = "none"  # Self negative loops are not shown
            else:
                mec = "orange"
            ax.plot(v1[0], v1[1] - 2*sep, "o", mec=mec, mew=2, mfc="none", markersize=10)
            continue

        plotted_pairs.append((a, b))

        if influence_graph.has_edge(b, a):
            b_reactions = influence_graph.edges[b, a]["reactions"]
            b_sign = influence_graph.edges[b, a]["sign"]
        else:
            b_reactions = []
            b_sign = 0

        nedges = len(reactions) + len(b_reactions)
        v1 = layout[a]
        v2 = layout[b]

        patchA = node_patches[a]
        patchB = node_patches[b]

        maxk = nedges//2 - (1 - nedges%2)*0.5
        kiter = np.arange(-maxk, maxk + 1)

        for k, R in zip(kiter, reactions):
            e = MutableEdge(ax,
                       nedges=nedges,
                       k=k,
                       name=R,
                       sign=sign,
                       v1=v1,
                       v2=v2,
                       patchA=patchA,
                       patchB=patchB,
                       sep=sep
            )
            edges[(a, b)][R] = e

        for k, R in zip(kiter + len(reactions), reversed(b_reactions)):
            e = MutableEdge(ax,
                       nedges=nedges,
                       k=-k,
                       name=R,
                       sign=b_sign,
                       v1=v2,
                       v2=v1,
                       patchA=patchB,
                       patchB=patchA,
                       sep=sep
            )
            edges[(b, a)][R] = e

    return MutableInfluenceGraph(ax, edges, max_edge_length)
