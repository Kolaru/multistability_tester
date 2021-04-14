import re

import networkx as nx
import numpy as np

from collections import OrderedDict
from itertools import product, chain
from numpy.linalg import det

from utils import pairs


def test_multistability(reaction_data):
    species, reactions = parse_reactions(reaction_data)
    contribution_graph = construct_contribution_graph(species, reactions)
    GI = construct_influence_graph(contribution_graph)
    cycles_info = retrieve_cycles_info(GI)
    return test_hoopings(reactions, cycles_info), GI


"""
    parse_species(s)

Parse the species in `s` given as a string `"nS"` where `n` is either a number or omitted, and `S` is the name
of the specy. Return `None` if passed `"0"`.
"""
def parse_species(s):
    if s == "0":
        return None

    n, S = re.findall(r'^(\d*)(\w+)', s)[0]

    if len(n) == 0:
        n = 1
    else:
        n = int(n)

    return n, S

"""
    parse_reactions(chem_network)

Parse reactions from a string describing a chemical network. Return the list of
species involved and a dictionnary describing the chemical reactions.
"""
def split_reactions(chem_network):
    # Parse all reactions and put them into a dictionnary
    react_raw = [s.strip() for s in chem_network.split("\n") if len(s.replace(" ", "")) > 0]  # Get a list of the lines
    react_str = [s.replace(" ", "") for s in react_raw]  # Remove all spaces
    react_str = [s for s in react_str if len(s) > 0]  # Remove empty lines

    reaction_names = []
    reaction_data = []

    for k, (raw, react) in enumerate(zip(react_raw, react_str)):
        nr = react.split(":")

        if len(nr) > 2:
            raise ValueError("all reactions may have at most one ':' between the "
                             "name of the reaction and the reaction itself. This "
                             "is not the case in reaction '{}'.".format(raw))

        if len(nr) == 2:
            named = True
            R, react = nr
            if R in reaction_names:
                raise ValueError("all reaction names must be unique. '{}' "
                                 "appears multiple time.")

            reaction_names.append(R)
        else:
            named = False
            R = "R{}".format(k)

        ba = react.split("<->")

        if len(ba) == 2:
            reversible = True
        else:
            reversible = False
            ba = react.split("->")

            if len(ba) != 2:
                raise ValueError("all reactions should have exactly one '->' or "
                                 "'<->'. This is not the case in the reaction "
                                 "'{}'.".format(raw))

        reaction_data.append(dict(
            before=ba[0],
            after=ba[1],
            reversible=reversible,
            name=R,
            named=named
        ))

    return reaction_data


def parse_reactions(reaction_data):
    species = set()
    reactions = OrderedDict()

    for k, data in enumerate(reaction_data):
        reactants = parse_side(data["before"])
        products = parse_side(data["after"])
        named = data["named"]
        R = data["name"]
        reversible=data["reversible"]

        if reversible:
            if named:
                names = (R + "+", R + "-")
            else:
                names = (R, "R-{}".format(k))
        else:
            names = (R,)

        spec = set(reactants.keys())
        spec.update(products.keys())
        species.update(spec)

        balance = {s:(products.get(s, 0) - reactants.get(s, 0)) for s in spec}

        # Forward and backward reactions are splitted if the reaction is reversible
        if len(names) == 2:
            Rf, Rb = names

            reactions[Rf] = dict(reactants=reactants,
                                 products=products,
                                 balance=balance)

            reactions[Rb] = dict(reactants=products,
                                 products=reactants,
                                 balance={s:-val for s, val in balance.items()})
        else:
            R = names[0]
            reactions[R] = dict(reactants=reactants,
                                products=products,
                                balance=balance)

    return sorted(list(species)), reactions


def parse_side(side):
    strings = side.split("+")
    complexes = [parse_species(s) for s in strings if s != "0"]
    return {c[1]:c[0] for c in complexes}

def construct_contribution_graph(species, reactions):
    G = nx.DiGraph()
    G.add_nodes_from(species)

    for R, reaction in reactions.items():
        for r in reaction["reactants"]:
            for b, val in reaction["balance"].items():
                if val != 0:
                    if G.has_edge(r, b):
                        G.edges[r, b]["contributions"][R] = val
                    else:
                        G.add_edge(r, b, contributions={R: val})

    return G


def construct_influence_graph(contribution_graph):
    GI = nx.DiGraph()  # Interaction graph
    GI.add_nodes_from(contribution_graph.nodes)

    for s1, s2, contributions in contribution_graph.edges(data="contributions"):
        vals = contributions.values()
        if all([v > 0 for v in vals]):
            sign = +1
        elif all([v < 0 for v in vals]):
            sign = -1
        else:
            sign = 0  # Impossible to determine the sign without knowing the kinetic

        GI.add_edge(s1, s2, sign=sign, reactions=list(contributions.keys()))

    return GI


def retrieve_cycles_info(GI):
    cycles = [tuple(c) for c in nx.simple_cycles(GI)] # Convert cycle to tuple to be able to use them as key
    cycles_info = []

    # Cycles are found as sequence of nodes, all possible edge combination
    # must be found for each cycle. The sign of each cycle do not depend on
    # the edges however.
    for cycle in cycles:
        paths = [[]]
        sign = 1
        for p in pairs(cycle):
            for k, path in enumerate(paths):
                # Replace each path/sign by a list of possible path/sign
                paths[k] = [path + [R] for R in GI.edges[p]["reactions"]]

            sign *= GI.edges[p]["sign"]

            # Flatten the lists
            paths = list(chain.from_iterable(paths))

        cycles_info.append(dict(cycle=cycle, paths=paths, sign=sign))

    return cycles_info


def extend_hooping(hooping, cycles):
    used = set(chain.from_iterable(hooping))
    compat = [(k, c) for k, c in enumerate(cycles) if len(used.intersection(c)) == 0]
    return [hooping + [c] for c in compat]


def test_hoopings(reactions, cycles):
    multistability = False

    cycles = sorted(cycles, key=lambda c: (c["sign"], len(c["cycle"])), reverse=True)

    n = 0
    for k, c in enumerate(cycles):
        # Since cycles are sorted, once we reach a cycle of negative sign, only negative sign cycles remain in the list.
        # Therefore, the hoopings build from them will never contain a positive sign cycle.
        if c["sign"] < 0:
            break

        queue = [[c]]

        while len(queue) > 0:
            hooping = queue.pop(0)
            species = list(chain.from_iterable([subcycle["cycle"] for subcycle in hooping]))

            # Find all possible combination of reactions
            for subpaths in product(*[subcycle["paths"] for subcycle in hooping]):
                n += 1

                # PERF using a numpy array of fixed size may help
                Rs = list(chain.from_iterable(subpaths))
                stoch = []
                for R in Rs:
                    stoch.append([reactions[R]["balance"].get(spec, 0) for spec in species])

                # TOFIX det always return a float, risk of imprecision breaking the code
                d = det(stoch)
                if d != 0:
                    multistability = True
                    return dict(
                        possible_multistability=True,
                        hooping=tuple(subcycle["cycle"] for subcycle in hooping),
                        path=subpaths,
                        det=d,
                        hoopings_tested=n,
                    )

            # PERF using a view on cycles may increase performance here
            queue.extend(extend_hooping(hooping, cycles[k+1:]))

    if not multistability:
        return dict(
            possible_multistability=False,
            hooping=None,
            path=None,
            det=None,
            hoopings_tested=n
        )

    return True
