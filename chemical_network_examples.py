import webbrowser

"""
    open_reference(network_example)

Open the reference for a given example in the default web browser.
"""
def open_reference(network_example):
    url = network_example["reflink"]

    if url is None:
        print("No link to reference available.")
    else:
        webbrowser.open(url)


double_negative_feedback = dict(
    network="""
        degA: A -> 0
        bindA: 2 A + OB <-> OB2A
        prodA: OA -> OA + A

        degB: B -> 0
        bindB: 2 B + OA <-> OA2B
        prodB: OB -> OB + B
    """,
    multistability=True,
    name="Double negative feedback loop",
    ref="Richard, Richard et al.",  # TODO find actual ref.
    reflink=None)

# Soliman 2013
ref_soliman2013 = ("Soliman, S. (2013). A stronger necessary condition for the "
                   "multistationarity of chemical reaction networks. Bulletin of "
                   "mathematical biology, 75(11), 2289-2303.")
reflink_soliman2013 = r"https://link.springer.com/article/10.1007/s11538-013-9893-7"

simple_enzyme_kinetics = dict(
    network="""
        E + S <-> EP
        EP -> E + P
    """,
    multistability=False,
    name="Simple enzyme kinetics",
    ref=ref_soliman2013,
    reflink=reflink_soliman2013)

two_step_enzyme_kinetics = dict(
    network="""
        M + K <-> MKq
        MK -> K + Mp
        K + Mp <-> MpK
        MpK -> K + Mpp
    """,
    multistability=False,
    name="Two step enzyme kinetics",
    ref=ref_soliman2013,
    reflink=reflink_soliman2013)

# Others
simple_opposite_influence = dict(
    network="""
        A + B -> C
        A + C -> D
    """,
    multistability=False,
    name="Simple opposite influence",
    ref=None,
    reflink=None)

# Siegal-Gaskins et al. 2011
ref_siegal_gaskins2011 = ("Siegal-Gaskins, D., Mejia-Guerra, M. K.,"
    "Smith, G. D., & Grotewold, E. (2011). Emergence of switch-like"
    "behavior in a large family of simple biochemical networks. PLoS"
    "computational biology, 7(5), e1002039.")
reflink_siegal_gaskins2011 = "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002039"

base_network = """
    X1 -> X1 + P1
    X2 -> X2 + P2
    P1 -> 0
    P2 -> 0
"""

reactions = dict(
    a="X1 + P1 <-> X1P1",
    b="X1 + P2 <-> X1P2",
    c="X2 + P1 <-> X2P1",
    d="X2 + P2 <-> X2P2",
    e="X1P1 -> X1P1 + P1",
    f="X1P2 -> X1P2 + P1",
    g="X2P1 -> X2P1 + P2",
    h="X2P2 -> X2P2 + P2",
    i="P1 + P1 <-> P1P1",
    j="P1 + P2 <-> P1P2",
    k="P2 + P2 <-> P2P2",
    l="X1 + P1P1 <-> X1P1P1",
    m="X1 + P1P2 <-> X1P1P2",
    n="X1 + P2P2 <-> X1P2P2",
    o="X2 + P1P1 <-> X2P1P1",
    p="X2 + P1P2 <-> X2P1P2",
    q="X2 + P2P2 <-> X2P2P2",
    r="X1P1P1 -> X1P1P1 + P1",
    s="X1P1P2 -> X1P1P2 + P1",
    t="X1P2P2 -> X1P2P2 + P1",
    u="X2P1P1 -> X2P1P1 + P2",
    v="X2P1P2 -> X2P1P2 + P2",
    w="X2P2P2 -> X2P2P2 + P2"
)

# From Fig. 5
minimal_codes = ["kqw", "bcdh", "bfjpv", "ckn", "abejp", "jmpsv",
                "ikno", "jknptv", "aejknp", "jkmnps", "dhjknp"]

minimal_bistable_networks = [
    dict(
        network=base_network + "\n" + "\n".join(r + ":" + reactions[r] for r in code),
        multistability=True,
        name="Minimal bistable network ({})".format(code),
        ref=ref_siegal_gaskins2011,
        reflink=reflink_siegal_gaskins2011
    ) for code in minimal_codes
]


examples = [double_negative_feedback,
            simple_enzyme_kinetics,
            two_step_enzyme_kinetics] + minimal_bistable_networks
