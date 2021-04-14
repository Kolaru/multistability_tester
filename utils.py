import itertools as itools

"""
    pairs(iter)

Return all consecutive pairs (including the last-first pair) of the iteraterable `iter`.

Inspired by itertools recipes.
"""
def pairs(iter):
    a, b = itools.tee(iter)
    b = itools.cycle(b)
    next(b, None)
    return zip(a, b)
