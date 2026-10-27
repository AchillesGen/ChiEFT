"""UFO loading, Dirac algebra and the evaluator for UFO Lorentz structures.

Conventions: metric diag(1,-1,-1,-1), Dirac representation, ubar u = 2m,
eps^{0123} = -1.  UFO legs are all outgoing, so an incoming particle enters a
vertex under its antiparticle name carrying minus its momentum.
"""
import cmath
import importlib
import itertools
import math
import os
import re
import sys

import numpy as np

HBARC2 = 389.379  # GeV^2 ub


# ---------------------------------------------------------------- UFO and card
def load_ufo(path):
    """(parameters, vertices, lorentz) modules of the UFO at `path`."""
    path = os.path.abspath(path)
    parent, pkg = os.path.split(path)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    return tuple(importlib.import_module(f"{pkg}.{m}")
                 for m in ("parameters", "vertices", "lorentz"))


def read_param_card(path):
    """{(BLOCK, (codes...)): value} from an SLHA param card, DECAY included."""
    out, block = {}, None
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line:
            continue
        words = line.split()
        key = words[0].upper()
        if key == "BLOCK":
            block = words[1].upper()
        elif key == "DECAY":
            block = None
            out[("DECAY", (int(words[1]),))] = float(words[2])
        elif block is not None:
            out[(block, tuple(int(w) for w in words[:-1]))] = float(words[-1])
    return out


def parameter_namespace(P, card=None, overrides=None):
    """Every UFO parameter resolved to a number.

    External parameters take their value from `overrides` (by name), else from
    the param card `card`, else the UFO default; internal ones are then
    evaluated from them.
    """
    card = read_param_card(card) if card else {}
    overrides = overrides or {}
    ns = {"cmath": cmath, "math": math, "complex": complex}
    for p in P.all_parameters:
        if p.nature == "external":
            key = (p.lhablock.upper(), tuple(p.lhacode))
            ns[p.name] = overrides.get(p.name, card.get(key, p.value))
    todo = [p for p in P.all_parameters if p.nature != "external"]
    while todo:
        rest = []
        for p in todo:
            try:
                ns[p.name] = eval(str(p.value), ns)
            except NameError:
                rest.append(p)
        if len(rest) == len(todo):
            raise RuntimeError("unresolved UFO parameters: %s"
                               % [p.name for p in rest][:6])
        todo = rest
    return ns


def declared_form_factors(L):
    """Names of the form-factor functions the UFO Lorentz structures call."""
    names = set()
    for l in L.all_lorentz:
        names |= set(re.findall(r"\b(F[A-Za-z0-9]*)\(\s*-?\s*P\(", l.structure))
    return sorted(names)


def re_(x):
    """Real part of a UFO parameter (internal ones come through complex)."""
    return float(x.real) if isinstance(x, complex) else float(x)


# ---------------------------------------------------------------- Dirac algebra
METRIC = np.diag([1.0, -1.0, -1.0, -1.0]).astype(complex)
G = np.diag(METRIC).real

_s = [np.array([[0, 1], [1, 0]], dtype=complex),
      np.array([[0, -1j], [1j, 0]], dtype=complex),
      np.array([[1, 0], [0, -1]], dtype=complex)]
_I2, _Z2 = np.eye(2, dtype=complex), np.zeros((2, 2), dtype=complex)

GAMMA = np.zeros((4, 4, 4), dtype=complex)
GAMMA[0] = np.block([[_I2, _Z2], [_Z2, -_I2]])
for _i in range(3):
    GAMMA[_i + 1] = np.block([[_Z2, _s[_i]], [-_s[_i], _Z2]])
GAMMA5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
ID4 = np.eye(4, dtype=complex)

SIGMA = np.zeros((4, 4, 4, 4), dtype=complex)        # (i/2)[gamma^mu, gamma^nu]
for _m in range(4):
    for _n in range(4):
        SIGMA[_m, _n] = 0.5j * (GAMMA[_m] @ GAMMA[_n] - GAMMA[_n] @ GAMMA[_m])

EPS = np.zeros((4, 4, 4, 4), dtype=complex)          # all indices up
for _perm in itertools.permutations(range(4)):
    EPS[_perm] = -np.linalg.det(np.eye(4)[list(_perm)])


def slash(v):
    return sum(v[m] * METRIC[m, m] * GAMMA[m] for m in range(4))


def dot(a, b):
    return complex(a @ METRIC @ b)


def uspinor(p, m, s):
    """Dirac spinor u(p, s), normalised to ubar u = 2m."""
    E = p[0].real
    chi = np.array([1.0, 0.0]) if s == 0 else np.array([0.0, 1.0])
    lower = (sum(p[i + 1] * _s[i] for i in range(3)) @ chi) / (E + m)
    return np.sqrt(E + m) * np.concatenate([chi, lower]).astype(complex)


def ubar(u):
    return np.conjugate(u) @ GAMMA[0]


# transverse photon polarisations, helicity -1 and +1, for q along +z
POL = [(-1.0 / np.sqrt(2)) * np.array([0, 1, 1j, 0], dtype=complex),
       (+1.0 / np.sqrt(2)) * np.array([0, 1, -1j, 0], dtype=complex)]


def kmom(s, m, MN):
    """Break-up momentum of sqrt(s) -> N(MN) + pi(m); 0 below threshold."""
    r = (s - MN ** 2 - m ** 2) ** 2 - 4.0 * m ** 2 * MN ** 2
    return 0.0 if (r <= 0.0 or s <= 0.0) else np.sqrt(r) / (2.0 * np.sqrt(s))


# ------------------------------------------------------------------ propagators
def prop_fermion(v, m):
    return 1j * (slash(v) + m * ID4) / (dot(v, v) - m ** 2)


def prop_scalar(v, m):
    return 1j / (dot(v, v) - m ** 2)


def prop_vector(v, M, W):
    """Unitary gauge, -i (g - v v / M^2) / (v^2 - M^2 + i M W); Feynman for M = 0."""
    den = dot(v, v) - M ** 2 + 1j * M * W
    D = np.zeros((4, 4), dtype=complex)
    for a in range(4):
        for b in range(4):
            D[a, b] = -1j * (METRIC[a, b] - (v[a] * v[b] / M ** 2 if M else 0.0)) / den
    return D


def rs_numerator(v, M):
    """Rarita-Schwinger numerator N[nu, eta] (4x4 blocks) of the spin-3/2 propagator."""
    sl = slash(v) + M * ID4
    out = np.zeros((4, 4, 4, 4), dtype=complex)
    for nu in range(4):
        for et in range(4):
            out[nu, et] = sl @ (-METRIC[nu, et] * ID4
                                + (1.0 / 3.0) * GAMMA[nu] @ GAMMA[et]
                                + (2.0 / (3.0 * M ** 2)) * v[nu] * v[et] * ID4
                                - (1.0 / (3.0 * M)) * (v[nu] * GAMMA[et] - GAMMA[nu] * v[et]))
    return out


def delta_line(V_left, V_right, v, M, D):
    """sum_{nu,eta} V_left[nu] (i N[nu,eta] / D) V_right[eta, mu], for every mu.

    V_left[nu] is the vertex on the left of the Delta line, V_right[eta, mu] the
    one on the right with the external boson index mu; either may carry mu."""
    N = rs_numerator(v, M)
    out = np.zeros((4, 4, 4), dtype=complex)
    for m in range(4):
        for nu in range(4):
            for et in range(4):
                L = V_left[nu, m] if V_left.ndim == 4 else V_left[nu]
                R = V_right[et, m] if V_right.ndim == 4 else V_right[et]
                out[m] += G[nu] * G[et] * L @ (1j * N[nu, et] / D) @ R
    return out


# ------------------------------------------------- Lorentz-structure evaluator
class T:
    """Tensor with labelled axes: ('L', k) Lorentz, ('S', k) spinor.

    Repeated labels are contracted (Lorentz with the metric), following the UFO
    convention that negative labels are summed over."""

    __slots__ = ("a", "axes")

    def __init__(self, a, axes):
        self.a = np.asarray(a, dtype=complex)
        self.axes = list(axes)

    def __mul__(self, o):
        if not isinstance(o, T):
            return T(self.a * o, self.axes)
        return T(np.tensordot(self.a, o.a, axes=0), self.axes + o.axes)._contract()

    __rmul__ = __mul__

    def __truediv__(self, o):
        return T(self.a / o, self.axes)

    def __neg__(self):
        return T(-self.a, self.axes)

    def __pow__(self, n):
        r = self
        for _ in range(int(n) - 1):
            r = r * self
        return r

    def __add__(self, o):
        if self.axes == o.axes:
            return T(self.a + o.a, self.axes)
        return T(self.a + np.transpose(o.a, [o.axes.index(x) for x in self.axes]),
                 self.axes)

    def __sub__(self, o):
        return self + (-o)

    def _contract(self):
        t = self
        while True:
            dup = next((x for x in t.axes if t.axes.count(x) == 2), None)
            if dup is None:
                return complex(t.a) if not t.axes else t
            i = t.axes.index(dup)
            j = t.axes.index(dup, i + 1)
            if dup[0] == "L":
                a = np.tensordot(t.a, METRIC, axes=([i], [0]))
                a = np.trace(a, axis1=j - 1 if j > i else j, axis2=a.ndim - 1)
            else:
                a = np.trace(t.a, axis1=i, axis2=j)
            axes = [x for n, x in enumerate(t.axes) if n not in (i, j)]
            if not axes:
                return complex(a)
            t = T(a, axes)


def eval_structure(struct, mom, ff):
    """Evaluate a UFO Lorentz structure; mom[i] is leg i's outgoing momentum."""
    ns = dict(
        P=lambda mu, n: T(mom[n], [("L", mu)]),
        Gamma=lambda mu, i, j: T(GAMMA, [("L", mu), ("S", i), ("S", j)]),
        Gamma5=lambda i, j: T(GAMMA5, [("S", i), ("S", j)]),
        Identity=lambda i, j: T(ID4, [("S", i), ("S", j)]),
        Metric=lambda mu, nu: T(METRIC, [("L", mu), ("L", nu)]),
        Sigma=lambda mu, nu, i, j: T(SIGMA, [("L", mu), ("L", nu), ("S", i), ("S", j)]),
        Epsilon=lambda a, b, c, d: T(EPS, [("L", a), ("L", b), ("L", c), ("L", d)]),
        ProjM=lambda i, j: T(0.5 * (ID4 - GAMMA5), [("S", i), ("S", j)]),
        ProjP=lambda i, j: T(0.5 * (ID4 + GAMMA5), [("S", i), ("S", j)]),
        complex=complex, cmath=cmath)
    ns.update(ff)
    return eval(struct, ns)


def vertex(V, names):
    """The unique UFO vertex with exactly these particles."""
    want = sorted(names)
    found = [v for v in V.all_vertices if sorted(p.name for p in v.particles) == want]
    if len(found) != 1:
        raise KeyError("%d UFO vertices with particles %s" % (len(found), names))
    return found[0]


def _evaluate(v, momout, ns, ff):
    names = [p.name for p in v.particles]
    legs = {i + 1: momout[n] for i, n in enumerate(names)}
    tot = None
    for (_, li), c in v.couplings.items():
        piece = eval_structure(v.lorentz[li].structure, legs, ff) * eval(str(c.value), dict(ns, **ff))
        tot = piece if tot is None else tot + piece
    return tot, names


def build(V, momout, ns, ff):
    """Vertex with fermions, as an array [RS Lorentz, vector Lorentz, row, col].

    momout maps each leg's particle name to its outgoing momentum.  Absent axes
    are omitted.  The row (ubar) index is the fermion leg with the larger PDG
    code, the column the other one."""
    v = vertex(V, list(momout))
    tot, names = _evaluate(v, momout, ns, ff)
    parts = v.particles
    ferm = [(i + 1, p) for i, p in enumerate(parts) if p.spin in (2, 4)]
    vec = [i + 1 for i, p in enumerate(parts) if p.spin == 3]
    rs = [i + 1 for i, p in enumerate(parts) if p.spin == 4]
    row = max(ferm, key=lambda x: x[1].pdg_code)[0]
    col = min(ferm, key=lambda x: x[1].pdg_code)[0]
    order = [("L", n) for n in rs] + [("L", n) for n in vec] + [("S", row), ("S", col)]
    return np.transpose(tot.a, [tot.axes.index(x) for x in order])


def build_bosonic(V, momout, ns, ff, vectors):
    """Vertex without fermions, as an array over the named vector legs, in order."""
    tot, names = _evaluate(vertex(V, list(momout)), momout, ns, ff)
    return np.transpose(tot.a, [tot.axes.index(("L", names.index(n) + 1)) for n in vectors])
