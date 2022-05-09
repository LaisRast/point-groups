#!/usr/bin/env python 
# coding: utf-8

## Generates decompositions of polar orbit polytopes of cyclic-type tubical
## groups $G$ when the image of the starting point is a rotation center of $G^h$.
## Outputs are json files which are stored in `data/`.
## These files are used in the online gallery.
## RUNTIME: 53680.75014424324s
##
## The cells are constructed using inexact calculations. This is fine for
## almost all cases, except some which we list in `exact_cases`. For these
## exceptional cases, we switch to exact calculations if running time is not
## so long. Otherwise we keep inexact calculations.
##
## NOTE: This script requires `PyNormaliz`. To install it, type
## `sage -i pynormaliz` in the terminal.

load("main_functions.sage")
from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
import json
import math

###############
## functions ##
###############
class AffineHullProjection:
    """
    Projection onto the affine hull of `F` which respects the "view from
    outside" convension.

    INPUT:
    - `F` -- Polyhedron: A 3-polytope in R^4, which is assumed to be the facet
      of a 4-polytope that contains the origin in its interior
    - `z_axis` -- tuple of two 4-lists (default: `None`): Endpoints of a line
      segment perpendicular to the affine hull of `F`.
    - `linear` -- boolean (default: `True`).

    OUTPUT:
    A function R^4 -> R^3: v |-> A(v) + b, which represents the projection
    onto the affine hull of `F`. If `linear`, then b is zero.
    The rows of the 3x4 orthonormal matrix A with the outer normal of `F`
    form a positive orthonormal basis.
    A and b are chosen such that the line segment defined by the pair of
    points in `z_axis` is mapped to a line parallel to the z-axis in R^3.
    """

    def __init__(self, F, z_axis=None, linear=True):
        ## Find an orthonormal basis
        affine_basis = F.an_affine_basis()
        M = matrix(AA, [v.vector() - affine_basis[0].vector() for v in affine_basis[1:]])
        A = M.gram_schmidt(orthonormal=True)[0]
        b = -A*vector(affine_basis[0])

        ## Ensure the "view from outside" convension
        # compute x4, the outer normal of F
        H_b = F.equations_list()[0][0]
        H_a = -vector(F.equations_list()[0][1:])
        x4 = H_a/H_a.norm()
        if not (H_a*vector([0,0,0,0]) <= H_b): x4 = -x4
        # Orient the basis using our convension
        if matrix([A[0], A[1], A[2], x4]).n().change_ring(QQ).det() > 0:
            A = matrix([A[0], A[1], A[2]])
        else:
            A = matrix([-A[0], A[1], A[2]])
        # Avoid nonrational numbers
        A = A.n().change_ring(QQ)

        ## Make sure that A maps the line in z_axis to the z-axis
        if z_axis:
            proj = lambda x: A*vector(x)+b
            v3 = proj(z_axis[0]) - proj(z_axis[1])
            v3 = v3/v3.norm()
            v2 = v3.cross_product(vector([1,2,3]))
            assert not v2.is_zero()
            v2 = v2/v2.norm()
            v1 = v2.cross_product(v3)
            v1 = v1/v1.norm()
            T = matrix([v1, v2, v3]).n().change_ring(QQ)
            assert T.det() > 0
            A = T*A
            b = T*b
        self.A = A.n().change_ring(QQ)
        self.b = vector([0,0,0]) if linear else b.n().change_ring(QQ)
    
    def __call__(self, x):
        return self.A*vector(x) + self.b

class CentralProjection:
    """
    Central projection on `F` which respects the "view from outside"
    convention. By default, this is inexact.

    INPUT:
    - `F` -- Polyhedron: A 3-polytope in R^4, which is assumed to be the facet
      of a 4-polytope that contains the origin in its interior
    - `projection_center` -- a 4-tuple.
    - `z_axis` -- tuple of two 4-lists (default: `None`): Endpoints of a line
      segment perpendicular to the affine hull of `F`.
    - `linear` -- boolean (default: `True`).

    OUTPUT:
    A function representing the central projection from `projection_center`
    onto `F`.

    REMARK:
    Adapted from ProjectionFuncSchlegel of sage >= 9.3.
    """
    def __init__(self, facet, projection_center, z_axis=None, linear=True):
        ineq = facet.equations()[0]
        self.full_A = ineq.A()
        self.full_b = ineq.b()
        self.projection_center = vector(RDF, projection_center)
        _proj = AffineHullProjection(facet, z_axis, linear=False)
        A, b = _proj.A, _proj.b
        A = matrix([-A[0], A[1], A[2]])
        if linear:
            self.proj = lambda x: A*vector(x) #+ b
        else:
            self.proj = lambda x: A*vector(x) + b

    def __call__(self, x):
        # The intersection of the segment with the facet
        # See Ziegler's "Lectures on Polytopes" p.133
        vx = vector(x)
        z = (self.full_b)
        a = -(self.full_A)
        y = self.projection_center
        preimage = y + ((z-a*y)/(a*vx-a*y))*(vx - y)
        return self.proj(preimage)

def facetope(v, gr, inexact=False):
    """
    The facet of the polar orbit polytope of P(v,gr) which correspond to v.
    The polar orbit polytope is scaled by ||v||**2.
    """
    scale = vector(v).norm()**2
    O = [vector(g.act(v)) for g in gr]
    v = vector(v)
    ieqs = [[1*scale] + list(-u) for u in O if (u.n() - v.n()).norm() > ALMOST_ZERO];
    eqns = [[1*scale] + list(-v)];

    if inexact:
        ieqs=[vector(u).n().change_ring(QQ) for u in ieqs]
        eqns=[vector(u).n().change_ring(QQ) for u in eqns]
        return Polyhedron(ieqs=ieqs, eqns=eqns, backend="normaliz")
    return Polyhedron(ieqs=ieqs, eqns=eqns, backend="normaliz")

def compact_json(text):
    """
    Remove spaces, keep only 5 digits after the comma and insert line breaks.
    eg: 12.1234567890 becomes 12.12345
    """
    line = ""
    lines = []
    IN_STRING = False
    num_digits = -1 # no fractional point has been seen
    for c in text:
        if c in " \n" and not IN_STRING:
            if len(line)>80:
                lines.append(line+"\n")
                line = ""
            # otherwise ignore
            num_digits = -1
            continue
        elif "0"<=c<="9" and not IN_STRING:
            if num_digits>5:
                continue # ignore extra digits
            if num_digits>=0:
                num_digits+=1
        else:
            num_digits = -1
            
        if c==".":
            num_digits = 0
        elif c=='"': # hope there is no '\"' or such troubles
            IN_STRING = not IN_STRING
        line += c
    if line:
        lines.append(line+"\n")
    return "".join(lines)
    
def identify_tubes(Cs, v, P):
    """
    Return which tube of O(v, Cs) is mapped to which vertex of P.
    """
    tubes_to_3d = {}
    Us = [vector(u).n() for u in P.vertices_list()]
    Os = [hopf_map(Cs[ii][0].act(v)) for ii in range(len(Cs))]

    for ii in range(len(Cs)):
        o = hopf_map(Cs[ii][0].act(v))
        for jj in range(len(Us)):
            u = Us[jj]
            if (o-u).norm() <= ALMOST_ZERO:
                tubes_to_3d[ii] = jj
    return tubes_to_3d

def angle_to_neighbor(m, s, t, n, ncells):
    """
    g = [exp(pi/m*p), exp(s*pi/t*i)]
    h = [1, exp(pi/n*i)]
    H = <g, h>
    """
    denom = lcm(m*t, n)
    gcd_, a, b = xgcd(denom*(-1/m+s/t), denom/n)

    assert a*(-1/m+s/t) + b*(1/n) == 2/ncells
    # assert gcd_/denom == 2/ncells 

    # return a*(1/m+s/t) + b(1/n)
    return 2*a/m, 2/ncells

def mod_2pi(angle):
    """
    Return `angle` mod 2*pi.
    """
    while angle.n() >= 2*pi.n():
        angle -= 2*pi
    while angle.n() < 0:
        angle += 2*pi
    return angle

def mod_2(number):
    """
    Return `number` mod 2*pi.
    """
    while number >= 2:
        number -= 2
    while number < 0:
        number += 2
    return number

def rangle_on_Cp(g, p):
    """
    Return the angle of the action of g on C_p as a right rotation.
    """
    langle = atan2(((g.l-g.l[0])/p)[0], g.l[0])
    # gr = g.r.to_exact_quaternion()
    # ranglex = atan2(((gr-gr[0])/i)[0], gr[0])
    rangle = g.r.alpha*pi if g.r.p == 0 else 0
    return mod_2pi(-langle + rangle)

def rotation_angle(g):
    """
    Return the angle of g - approximated
    """
    langle = atan2(((g.l-g.l[0])/p)[0], g.l[0])
    # gr = g.r.to_exact_quaternion()
    # ranglex = atan2(((gr-gr[0])/i)[0], gr[0])
    rangle = g.r.alpha*pi if g.r.p == 0 else 0
    sangle1 = mod_2pi(langle+rangle)
    sangle2 = nearby_rational_angle(sangle1)*pi
    assert sangle1.n() - sangle2.n() <= ALMOST_ZERO
    return sangle2/pi

def number_to_latex(x):
    """
    x is rational.
    """
    if x.denominator() == 1: return f"{x}"
    if x.numerator() < 0:
        return f"-\\frac{{{-x.numerator()}}}{{{x.denominator()}}}"
    return f"\\frac{{{x.numerator()}}}{{{x.denominator()}}}"

## make FDR(a,b) and FDR(-a,-b) the same
def new_eq(self, other):
    lminus1 = quat_minus1 if hasattr(self.l, "alpha") else -1
    rminus1 = quat_minus1 if hasattr(self.r, "alpha") else -1
    return (self.l == other.l and self.r == other.r and self.star == other.star) or \
           (self.l == lminus1*other.l and self.r == rminus1*other.r and self.star == other.star)
def new_hash(self):
    lminus1 = quat_minus1 if hasattr(self.l, "alpha") else -1
    rminus1 = quat_minus1 if hasattr(self.r, "alpha") else -1
    if self.l < lminus1*self.l:
        return hash((lminus1*self.l,rminus1*self.r,self.star))
    return hash((self.l,self.r,self.star))
FDR.__eq__ = new_eq
FDR.__hash__ = new_hash

##########
## Hopf ##
##########
class HopfBundle:
    """
    The Hopf bundle $H_{q_0}$ or $H_^{q_0}$.
    """
    def __init__(self, q0, left=True):
        self.q0 = q0
        self.left = left

    def map(self, x):
        if self.left:
            return x * self.q0 * x.conjugate()
        return x.conjugate() * self.q0 * x

hopf_map = lambda x: vector(HopfBundle(i, left=True).map(Q(x)).coefficient_tuple()[1:])

class GreatCircle:
    """
    The great circle $K_p^q$.
    """
    def __init__(self, p, q):
        self.p = p
        self.q = q
        self.simple_rotation = FDR(p, q)

    def an_element(self, approx=False):
        """
        Return the point rot in $K_p^q$
        such that [rot] is the rotation with the smallest angle.
        """
        if self.p == self.q: return Q(1)
        if self.p == -self.q: return Q(-1)
        # formula: q*p = -dot(q, p) + cross(p, q)
        axis = self.q*self.p - (self.q*self.p)[0]
        axis = axis/vector(axis).norm()
        cos_angle = -(self.q*self.p)[0]
        if approx: cos_angle = QQ(cos_angle.n())
        sin_half_angle = sqrt((1-cos_angle)/2)
        cos_half_angle = sqrt((1+cos_angle)/2)
        if approx:
            sin_half_angle = QQ(sqrt((1-cos_angle)/2).n())
            cos_half_angle = QQ(sqrt((1+cos_angle)/2).n())
        rot = cos_half_angle + sin_half_angle*axis
        if not approx:
            assert self.contains(rot)
        return rot

    def contains(self, rot):
        """
        Test if `rot` lies on $K_p^q$.
        """
        return FDR(rot, rot).act(self.p) == self.q

    def plot(self, proj=None, n_points=100, **kwd):
        """
        Plot the circle into 3d using the function $proj$.
        """
        # if proj is None:
        #     proj = lambda x: x
        circle_points = []
        an_element_inexact = Q_in(list(self.an_element(approx=True)))
        q_inexact = Q_in(list(self.q))
        for ii in range(n_points):
            point = an_element_inexact*(cos(2*ii*pi.n()/n_points) + \
                                        sin(2*ii*pi.n()/n_points)*q_inexact)
            circle_points.append(proj(point))
        return line3d(circle_points + [circle_points[0]], **kwd)

###############
## constants ##
###############
phi = (1+AA(sqrt(5)))/2
Phi = 1/phi
DISCONNECTED = "the cells of a tube are disconnected from each other.\n"

############
## inputs ##
############
verbose = True
n_figures = 15 # multiple of 3 !

#####################################
## 3-dimensional polyhedral groups ##
#####################################
grp_2I = Q_group([i_I, w]).q
grp_2O = Q_group([i_O, w]).q
grp_2T = Q_group([i, w]).q
grp_pI = [FDR(g, g) for g in grp_2I]
grp_pO = [FDR(g, g) for g in grp_2O]
grp_pT = [FDR(g, g) for g in grp_2T]
grp_pC = lambda n: [FDR(g, g) for g in grp_2C(n)]
grp_pD = lambda n: [FDR(g, g) for g in grp_2D(n)]
grp_pmI = [FDR(g, g) for g in grp_2I] + [FDR(-g, g) for g in grp_2I]
grp_pmO = [FDR(g, g) for g in grp_2O]+ [FDR(-g, g) for g in grp_2O]
grp_pmT = [FDR(g, g) for g in grp_2T]+ [FDR(-g, g) for g in grp_2T]

#################
## groups data ##
#################
grps = [
    {
    "group": lambda n: mk_group([FDR(i_I, quat_1), FDR(w, quat_1), FDR(1, quat_2D(1/n))]),
    "gr_name": lambda n: f"IxC{n}" if n > 0 else "IxCn",
    "gr_latex": lambda n: f"$\\pm[I\\times C_{n}]$" if n > 0 else "$\\pm[I\\times C_n]$",
    "grh": grp_pI,
    "grh_latex": "${+I}$",
    "centers":
        [
            {"type": "5-fold",
             "ntubes": "12",
             "point": (0, Phi, 1),
             "pt_latex": "$\\frac1{\\sqrt{\\varphi^2+1}}(0, 1, \\varphi)$, where $\\varphi=\\frac{1+\\sqrt5}{2}$",
             "ncells": lambda n: lcm(2*n, 10),
             "ncells_latex": "$\\mathrm{lcm}(2n, 10)$",
             "tubes_to_3d": {0: 7, 1: 11, 2: 9, 3: 5, 4: 10, 5: 1, 6: 3, 7: 8, 8: 6, 9: 2, 10: 4, 11: 0},
             "element_g": -w*i_I,
             "H_gens": lambda n: [FDR(-w*i_I, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega i_I, 1], [1, e_n] \\rangle$",
             "alternate_grs": r"$\pm[I\times D_{2n}]$",
             "H_rotations": lambda n: (5, 0, 1, n),
            },
            {"type": "3-fold",
             "ntubes": "20",
             "point": (-1, -1, -1),
             "pt_latex": "$\\frac1{\\sqrt3}(-1, -1, -1)$",
             "ncells": lambda n: lcm(2*n, 6),
             "ncells_latex": "$\\mathrm{lcm}(2n, 6)$",
             "tubes_to_3d": {0: 2, 1: 6, 2: 0, 3: 12, 4: 8, 5: 1, 6: 4, 7: 14, 8: 9, 9: 3, 10: 10, 11: 16, 12: 7, 13: 5, 14: 19, 15: 18, 16: 15, 17: 11, 18: 13, 19: 17},
             "element_g": -w,
             "H_gens": lambda n: [FDR(-w, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega, 1], [1, e_n] \\rangle$",
             "alternate_grs": "$\\pm[I\\times D_{2n}]$",
             "H_rotations": lambda n: (3, 0, 1, n),
            },
            {"type": "2-fold",
             "ntubes": "30",
             "point": (1, Phi, Phi+1),
             "pt_latex": "$\\frac12(1, \\frac1\\varphi, \\varphi)$, where $\\varphi=\\frac{1+\\sqrt5}{2}$",
             "ncells": lambda n: lcm(2*n, 4),
             "ncells_latex": "$\\mathrm{lcm}(2n, 4)$",
             "tubes_to_3d": {0: 24, 1: 28, 2: 15, 3: 20, 4: 29, 5: 22, 6: 8, 7: 16, 8: 27, 9: 26, 10: 18, 11: 6, 12: 12, 13: 19, 14: 23, 15: 25, 16: 10, 17: 2, 18: 4, 19: 11, 20: 21, 21: 17, 22: 13, 23: 0, 24: 3, 25: 7, 26: 14, 27: 9, 28: 1, 29: 5},
             "element_g": i_I,
             "H_gens": lambda n: [FDR(i_I, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [i_I, 1], [1, e_n] \\rangle$",
             "alternate_grs": "$\\pm[I\\times D_{2n}]$",
             "H_rotations": lambda n: (2, 0, 1, n),
             "polyhedron": "dodecahedron",
            },
        ]
    },
    {
    "group": lambda n: mk_group([FDR(i_O, quat_1), FDR(w, quat_1), FDR(1, quat_2D(1/n))]),
    "gr_name": lambda n: f"OxC{n}" if n > 0 else "OxCn",
    "gr_latex": lambda n: f"$\\pm[O\\times C_{n}]$" if n > 0 else "$\\pm[O\\times C_n]$",
    "grh": grp_pO,
    "grh_latex": "${+O}$",
    "centers":
        [
            {"type": "4-fold",
             "ntubes": "6",
             "point": (0, 1, 0),
             "pt_latex": "$(0, 1, 0)$",
             "ncells": lambda n: lcm(2*n, 8),
             "ncells_latex": "$\\mathrm{lcm}(2n, 8)$",
             "tubes_to_3d": {0: 4, 1: 3, 2: 5, 3: 0, 4: 1, 5: 2},
             "element_g": -w*i_O,
             "H_gens": lambda n: [FDR(-w*i_O, quat_1), FDR(1, quat_2D(1/n))],
             "alternate_grs": "$\\pm[O\\times D_{2n}]$",
             "H_latex": "$\\langle [-\\omega i_O, 1], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (4, 0, 1, n),
            },
            {"type": "3-fold",
             "ntubes": "8",
             "point": (-1, -1, -1),
             "pt_latex": "$\\frac1{\\sqrt3}(-1, -1, -1)$",
             "ncells": lambda n: lcm(2*n, 6),
             "ncells_latex": "$\\mathrm{lcm}(2n, 4)$",
             "tubes_to_3d": {0: 0, 1: 4, 2: 2, 3: 5, 4: 1, 5: 6, 6: 3, 7: 7},
             "element_g": -w,
             "H_gens": lambda n: [FDR(-w, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega, 1], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (3, 0, 1, n),
             "alternate_grs": "$\\pm[O\\times D_{2n}]$",
            },
            {"type": "2-fold",
             "ntubes": "12",
             "point": (0, 1, 1),
             "pt_latex": "$\\frac1{\\sqrt2}(0, 1, 1)$",
             "ncells": lambda n: lcm(2*n, 4),
             "ncells_latex": "$\\mathrm{lcm}(2n, 4)$",
             "tubes_to_3d": {0: 7, 1: 10, 2: 3, 3: 11, 4: 5, 5: 2, 6: 6, 7: 9, 8: 8, 9: 0, 10: 1, 11: 4},
             "H_rotations": lambda n: (2, 0, 1, n),
             "element_g": i_O,
             "H_gens": lambda n: [FDR(i_O, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [i_O, 1], [1, e_n] \\rangle$",
             "alternate_grs": "$\\pm[O\\times D_{2n}]$",
            },
        ]
    },
    {
    "group": lambda n: mk_group([FDR(i, quat_1), FDR(w, quat_1), FDR(1, quat_2D(1/n)), FDR(i_O, quat_2D(1/(2*n)))]),
    "gr_name": lambda n: f"hOxC{2*n}" if n > 0 else "hOxC2n",
    "gr_latex": lambda n: f"$\\pm\\frac12[O\\times C_{{{2*n}}}]$" if n > 0 else "$\\pm\\frac12[O\\times C_{2n}]$",
    "grh": grp_pO,
    "grh_latex": "${+O}$",
    "centers":
        [
            {"type": "4-fold",
             "ntubes": "6",
             "point": (0, 1, 0),
             "pt_latex": "$(0, 1, 0)$",
             "ncells": lambda n: 8*n/gcd(n-2, 4),
             "ncells_latex": "$\\frac{8n}{\\gcd(n-2, 4)}$",
             "tubes_to_3d": {0: 4, 1: 1, 2: 3, 3: 2, 4: 5, 5: 0},
             "element_g": -w*i_O,
             "H_gens": lambda n: [FDR(-w*i_O, quat_2D(1/(2*n))), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega i_O, e_{2n}], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (4, 1, 2*n, n),
             "alternate_grs": "$\\pm\\frac12[O\\times \\overline{O}_{4n}]$",
            },
            {"type": "3-fold",
             "ntubes": "8",
             "point": (-1, -1, -1),
             "pt_latex": "$\\frac1{\\sqrt3}(-1, -1, -1)$",
             "ncells": lambda n: lcm(2*n, 6),
             "ncells_latex": "$\\mathrm{lcm}(2n, 6)$",
             "tubes_to_3d": {0: 0, 1: 3, 2: 4, 3: 5, 4: 7, 5: 2, 6: 6, 7: 1},
             "H_rotations": lambda n: (3, 0, 1, n),
             "element_g": -w,
             "H_gens": lambda n: [FDR(-w, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega, 1], [1, e_n] \\rangle$",
             "alternate_grs": "$\\pm\\frac12[O\\times \\overline{D}_{4n}]$",
            },
            {"type": "2-fold",
             "ntubes": "12",
             "point": (0, 1, 1),
             "pt_latex": "$\\frac1{\\sqrt2}(0, 1, 1)$",
             "ncells": lambda n: lcm(2*n, 4*n/gcd(n-1, 8)),
             "ncells_latex": r"$\frac{4n}{\gcd(n-1, 2)}$",
             "tubes_to_3d": {0: 7, 1: 4, 2: 10, 3: 1, 4: 9, 5: 11, 6: 3, 7: 2, 8: 0, 9: 8, 10: 5, 11: 6},
             "element_g": i_O,
             "H_gens": lambda n: [FDR(i_O, quat_2D(1/(2*n))), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [i_O, e_{2n}], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (2, 1, 2*n, n),
             "alternate_grs": "$\\pm\\frac12[O\\times \\overline{D}_{4n}]$",
            },
        ]
    },
    {
    "group": lambda n: mk_group([FDR(i, quat_1), FDR(w, quat_1), FDR(1, quat_2D(1/n))]),
    "gr_name": lambda n: f"TxC{n}" if n > 0 else "TxCn",
    "gr_latex": lambda n: f"$\\pm[T\\times C_{n}]$" if n > 0 else "$\\pm[T\\times C_n]$",
    "grh": grp_pT,
    "grh_latex": "${+T}$",
    "centers":
        [
            {"type": "3-fold (type I)",
             "ntubes": "4",
             "point": (-1, -1, -1),
             "pt_latex": "$\\frac1{\\sqrt3}(-1, -1, -1)$",
             "ncells": lambda n: lcm(2*n, 6),
             "ncells_latex": "$\\mathrm{lcm}(2n, 6)$",
             "tubes_to_3d": {0: 0, 1: 1, 2: 2, 3: 3},
             "element_g": -w,
             "H_gens": lambda n: [FDR(-w, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega, 1], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (3, 0, 1, n),
             "alternate_grs": "$\\pm\\frac12[O\\times D_{2n}]$",
            },
            # {"type": "3-fold (type II)",
            #  "ntubes": "4",
            #  "point": (1, 1, 1),
            #  "pt_latex": "$\\frac1{\\sqrt3}(1, 1, 1)$",
            #  "ncells": lambda n: lcm(2*n, 6),
            #  "ncells_latex": "$\\mathrm{lcm}(2n, 6)$",
            #  "tubes_to_3d": {0: 3, 1: 2, 2: 1, 3: 0},
            #  "element_g": -w^2,
            #  "H_gens": lambda n: [FDR(-w^2, quat_1), FDR(1, quat_2D(1/n))],
            #  "H_latex": "$\\langle [-\\omega^2, 1], [1, e_n] \\rangle$",
            #  "H_rotations": lambda n: (3, 0, 1, n),
            #  "alternate_grs": "$\\pm\\frac12[O\\times D_{2n}]$",
            # },
            {"type": "2-fold",
             "ntubes": "6",
             "point": (1, 0, 0),
             "pt_latex": "$(1, 0, 0)$",
             "ncells": lambda n: lcm(2*n, 4),
             "ncells_latex": "$\\mathrm{lcm}(2n, 4)$",
             "tubes_to_3d": {0: 5, 1: 4, 2: 1, 3: 3, 4: 2, 5: 0},
             "element_g": i,
             "H_gens": lambda n: [FDR(i, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [i, 1], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (2, 0, 1, n),
             "alternate_grs": "$\\pm[T\\times D_{2n}]$ and $\\pm\\frac12[O\\times D_{2n}]$ (also their common supergroup $\\pm[O\\times D_{2n}]$)\nif $n \\equiv 0 \\mod 4$, else $\\pm[T\\times D_{2n}]$ (and its supergroup $\\pm\\frac12[O\\times\\overline{D}_{4n}]$)",
            },
        ]
    },
    {
    "group": lambda n: mk_group([FDR(i, quat_1), FDR(w, quat_2D(1/(3*n))), FDR(1, quat_2D(1/n))]),
    "gr_name": lambda n: f"tTxC{3*n}" if n > 0 else "tTxC3n",
    "gr_latex": lambda n: f"$\\pm\\frac13[T\\times C_{{{3*n}}}]$" if n > 0 else "$\\pm\\frac13[T\\times C_{3n}]$",
    "grh": grp_pT,
    "grh_latex": "${+T}$",
    "centers":
        [
            {"type": "3-fold (type I)",
             "ntubes": "4",
             "point": (-1, -1, -1),
             "pt_latex": "$\\frac1{\\sqrt3}(-1, -1, -1)$",
             # "ncells": lambda n: lcm(2*n, 6*n/gcd(n-1, 6)),
             "ncells": lambda n: 6*n/gcd(n-1, 3),
             # "ncells_latex": "$\\mathrm{lcm}(2n, \\frac{6n}{\\gcd(n-1, 6)})$",
             "ncells_latex": "$\\frac{6n}{\\gcd(n-1, 3)}$",
             "tubes_to_3d": {0: 0, 1: 1, 2: 2, 3: 3},
             "element_g": -w,
             "H_gens": lambda n: [FDR(-w, quat_2D(1/(3*n))), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega, e_{3n}], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (3, 1, 3*n, n),
             "alternate_grs": "$\\pm\\frac16[O\\times D_{6n}]$ (and its supergroup $\\pm\\frac12[O\\times D_{6n}]$ if $n\\not\\equiv 1\\mod 3$)",
            },
            {"type": "3-fold (type II)",
             "ntubes": "4",
             "point": (1, 1, 1),
             "pt_latex": "$\\frac1{\\sqrt3}(1, 1, 1)$",
             # "ncells": lambda n: lcm(2*n, 6*n/gcd(n-2, 12)),
             "ncells": lambda n: 6*n/gcd(n-2, 3),
             # "ncells_latex": "$\\mathrm{lcm}(2n, \\frac{6n}{\\gcd(n-2, 12)})$",
             "ncells_latex": "$\\frac{6n}{\\gcd(n-2, 3)}$",
             "tubes_to_3d": {0: 3, 1: 2, 2: 1, 3: 0},
             "H_rotations": lambda n: (3, 2, 3*n, n),
             "element_g": -w^2,
             "H_gens": lambda n: [FDR(-w^2, quat_2D(1/(3*n))^2), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [-\\omega^2, e_{3n}^2], [1, e_n] \\rangle$",
             "alternate_grs": "$\\pm\\frac16[O\\times D_{6n}]$ (and its supergroup $\\pm\\frac12[O\\times D_{6n}]$ if $n\\not\\equiv 2\\mod 3$)",
            },
            {"type": "2-fold",
             "ntubes": "6",
             "point": (1, 0, 0),
             "pt_latex": "$(1, 0, 0)$",
             "ncells": lambda n: lcm(2*n, 4),
             "ncells_latex": "$\\mathrm{lcm}(2n, 4)$",
             "tubes_to_3d": {0: 5, 1: 4, 2: 1, 3: 3, 4: 2, 5: 0},
             "element_g": i,
             "H_gens": lambda n: [FDR(i, quat_1), FDR(1, quat_2D(1/n))],
             "H_latex": "$\\langle [i, 1], [1, e_n] \\rangle$",
             "H_rotations": lambda n: (2, 0, 1, n),
             "alternate_grs": "$\\pm\\frac16[O\\times D_{6n}]$",
            },
        ]
    },
]

#######################
## Exceptional Cases ##
#######################
## Cases to be constructed using exact calculations.
## Commented ones produce wrong cells but they take a very long running time.
## Thus, we correct them by hand in the paper.
exact_cases = [
    "IxC75-5", # inexact group
    "IxC4-2",
    "OxC8-4", "OxC16-4", "OxC24-4",
    "OxC32-4", # 140.714296s
    "OxC40-4", # 281.241646s
    "OxC48-4", # 251.436913s
    "OxC18-3", # 922.978478s
    "OxC36-3", # 2905.743811s
    "OxC8-2",
    "hOxC4-4", "hOxC12-4", 
    "hOxC18-3", 
    # "hOxC54-3", # forever
    # "hOxC90-3", # forever
    "TxC36-3", # 1107.957481s; projected in 153.349783s
    "TxC18-3",
    "TxC18-3p", # 436.949342s
    # "TxC36-3p", # forever
    "tTxC3-3",
    "tTxC12-3", # 1064.774214s
    # "tTxC30-3", # 4977.961172s # sometimes does not work
    # "tTxC48-3", # forever
    # "tTxC66-3", # forever
    "tTxC6-3p",
    "tTxC24-3p",
    # "tTxC42-3p" # forever
    # "tTxC60-3p" # forever
]

## groups to be constructed using inexact FDR's
## since sage cannot compute some values, like e(75), exactly
inexact_groups = [
    "IxC75-5",
]


###############
## main loop ##
###############
start_overall = time()
for grp in grps:
    gr_func = grp["group"]
    gr_name_func = grp["gr_name"]
    gr_latex_func = grp["gr_latex"]
    gr_latex = gr_latex_func(0)
    centers = grp["centers"]
    grh = grp["grh"]
    grh_latex = grp["grh_latex"]

    ## center loop
    for center in centers[:]:
        pt = center["point"]
        pt_typ = center["type"]
        ntubes = center["ntubes"]
        pt_latex = center["pt_latex"]
        ncells_func = center["ncells"]
        ncells_latex = center["ncells_latex"]
        element_g = center["element_g"]
        H_rotations_func = center["H_rotations"]
        H_gens_func = center["H_gens"]
        H_latex = center["H_latex"]
        alternate_grs = center["alternate_grs"]
        
        k = pt_typ[0]
        if "type II" in pt_typ: k += "p"

        ## find starting point
        p = Q([0]+list(pt))/AA(vector(pt).norm())
        v = GreatCircle(p, i).an_element()
    
        if verbose:
            verbose_string = f"## BEGIN GROUP: {gr_name_func(0)}; CENTER: {pt_typ} ##"
            print("#" * len(verbose_string))
            print(verbose_string)
            print("#" * len(verbose_string))

        ## prepare gallery dict
        gallery = {"cells": []}
        gallery["information"] = {
            "gr_latex": f"{gr_latex}",
            "center_type": pt_typ,
            "k": k,
            "ntubes": int(ntubes),
        }
        
        ## sort ns according to ncells
        ncells_and_ns = {}
        for n in range(1, 1000):
            ncells=ncells_func(n)
            if ncells not in ncells_and_ns: ncells_and_ns[ncells] = []
            ncells_and_ns[ncells].append(n)
 
        ## sometimes there is not much space: decrease by 1 row
        n_figures_current = n_figures
        if gr_name_func(0)+"-"+k in ["hOxC2n-2","hOxC2n-4","tTxC3n-3","tTxC3n-3p","tTxC3n-2"]:
            n_figures_current -= 3 

        ## projection center should be computed once for each center 
        ## and it should be the same for all 
        compute_proj_center = True

        ## centers loop
        for ncells in sorted(list((ncells_and_ns.keys())))[:n_figures_current]:
            
            ## take the smallest n
            ns = sorted(ncells_and_ns[ncells])
            n = ns[0]

            ## name of the group
            gr_name = gr_name_func(n)
            if verbose: print(f"## BEGIN TUBE: group {gr_name}, {k}-fold center, n={n}")

            ## generate group
            start = time()
            gr = gr_func(n)
            runtime = f" constructed in {time() - start:f}s"
            if verbose: print(f"## GROUP: {gr_name_func(n)}"+runtime)

            ## angle to neighbor
            H = mk_group(H_gens_func(n))
            neighbor_gs = [g for g in H if abs(2*pi.n()/ncells - rangle_on_Cp(g, p).n()).n()<ALMOST_ZERO]
            to_neighbor_angle = angle_to_neighbor(*H_rotations_func(n), ncells)
            assert len(neighbor_gs) > 0
            assert mod_2(sum(to_neighbor_angle)) in [rotation_angle(g) for g in neighbor_gs]

            if len(neighbor_gs) == 1:
                to_neighbor_angle_latex = f"$({number_to_latex(mod_2(to_neighbor_angle[0])/2)} + "
                to_neighbor_angle_latex += f"{number_to_latex(to_neighbor_angle[1]/2)})\\cdot 2\\pi$"

            else:
                if to_neighbor_angle[0] == 0:
                    assert set([mod_2(2*k/H_rotations_func(n)[0] + to_neighbor_angle[1]) for k in range(10)]) == set([mod_2(rotation_angle(g)) for g in neighbor_gs])
                    to_neighbor_angle_latex = f"$(\\frac{{k}}{{{H_rotations_func(n)[0]}}} + "
                    to_neighbor_angle_latex += f"{number_to_latex(to_neighbor_angle[1]/2)})\\cdot 2\\pi$"
                else:
                    assert len(neighbor_gs) == 2 # only happens when gr=hOxC2n
                    assert set([mod_2((2*k+1)/2 + to_neighbor_angle[1]) for k in range(10)]) == set([mod_2(rotation_angle(g)) for g in neighbor_gs])
                    to_neighbor_angle_latex = f"$(\\frac{{2k+1}}{{4}} + "
                    to_neighbor_angle_latex += f"{number_to_latex(to_neighbor_angle[1]/2)})\\cdot 2\\pi$"
 

            ## construct facet
            inexact = True
            if gr_name_func(n) + "-" + k in exact_cases:
                inexact = False
            if verbose: print(f"CELL: inexact = {inexact}")

            if verbose: print(f"CELL: constructing...")
            start = time()
            if gr_name_func(n) + "-" + k in inexact_groups:
                v_in = Q_in(list(v))
                gr_in = [FDR(Q_in(list(g.l)), g.r.to_quaternion()) for g in gr]
                inexact = True
                if verbose: print(f"CELL: inexact = {inexact}")
                F = facetope(v_in, gr_in, inexact)
            else:
                F = facetope(v, gr, inexact)
            if verbose: print(f"CELL: constructed in {time() - start:f}s")

            ###############
            ## cell data ##
            ###############

            ## collect cell data for the gallery
            cell_data = {"points": []}
            cell_data["information"] = {
                "ns": ns,
                "ntubes": int(ntubes),
                "ncells_on_tube": int(ncells),
                "ncells_all": int(ntubes)*int(ncells),
                "inexact": inexact,
                "fname": f"{gr_name_func(0)}/{k}/{int(ncells)*int(ntubes)}cells_{ntubes}tubes"
            }
            
            ## projection
            z_axis = [v, FDR(1, e(14)).act(v)]
            proj = AffineHullProjection(F, z_axis)

            ## collect face
            F_faces = []
            for face_equation in F.inequality_generator():
                face_vertices = cyclic_sort_vertices_2d([vert for vert in face_equation.incident()])
                vert_0, vert_1, vert_2 = [proj(vector(v)) for v in face_vertices[:3]]
                face_normal = (vert_2 - vert_0).cross_product(vert_1 - vert_0)
                if face_normal.dot_product(vert_0) < 0:
                    face_vertices.reverse()
                F_faces.append([F.vertices().index(vert) for vert in face_vertices])
            cell_data["faces"] = F_faces
            F_vertices = [Q(p) for p in F.vertices_list()]

            ## collect vertices
            for pp in F_vertices:
                pp_proj  =tuple(proj(pp.coefficient_tuple()))
                cell_data["points"].append(tuple([float(x) for x in pp_proj]))
            
            gallery["cells"].append(cell_data)
            if verbose: print(f"CELL: cell data collected\n")

            ###########
            ## Tubes ##
            ###########
            start_constructing = time()

            ## subgroup and its cosets
            sgr = [g for g in gr if g.l in mk_group([element_g])]
            assert len(gr)/len(sgr) == int(ntubes)
            Cs = cosets(gr, sgr)

            ## remove duplications:
            ## elements of the subgroup which map v to the same point
            if len(sgr) != ncells:
                if gr_name_func(0) in ['OxCn', 'hOxC2n'] and ncells == 4*n:
                    new_sgr = []
                    visited_on_Cp = []
                    CS0_indices = []
                    for ii in range(len(Cs[0])):
                        g = Cs[0][ii]
                        gv = g.act(v)
                        if gv in visited_on_Cp: continue
                        visited_on_Cp.append(gv)
                        CS0_indices.append(ii)
                else:
                    CS0_indices = [ii for ii in range(len(Cs[0])) if Cs[0][ii].l in [-1, 1]]
                Cs = [[C[ii] for ii in CS0_indices] for C in Cs]
            assert len(Cs[0]) == ncells

            ## which tubes are mapped to which points in S^2
            P = Polyhedron([g.act(p).coefficient_tuple()[1:] for g in grh], backend="normaliz")
            assert identify_tubes(Cs, v, P) == center["tubes_to_3d"]

            ## projection
            z_axis = [v, FDR(1, e(14)).act(v)]
            if compute_proj_center:
                proj_center = 1.05*F.center()
                compute_proj_center = False
            proj_sch = CentralProjection(F, proj_center, z_axis=z_axis)

            ## sort point depending on the distance in R4 from -v.
            v_vector = vector(v).n()
            ## sort Cs
            for ii in range(len(Cs)):
                Cs[ii] = sorted(Cs[ii], key=lambda g: (vector(g.act(v)).n()+v_vector).norm())
            ## obtain cells order
            all_vertices = []
            for C in Cs:
                for g in C:
                    all_vertices.append(vector(g.act(v)).n())
            cells_order = sorted([ii for ii in range(len(all_vertices))], key=lambda ii: (all_vertices[ii]+v_vector).norm())
            cells_order = [[ii//int(ncells), ii%int(ncells)] for ii in cells_order]
    
            ## tubes data
            data =  {"points": [], "tubes": {}}
            data["information"] = {
                "gr_latex": f"{gr_latex}",
                "gr_latex_ns": [f"{gr_latex_func(n)}" for n in ns],
                "ns": ns,
                "to_neighbor_angle_latex": to_neighbor_angle_latex,
                "ntubes": int(ntubes),
                "ncells_on_tube": int(ncells),
                "ncells_all": int(ntubes)*int(ncells),
                "center_type": pt_typ,
                "k": k,
                "inexact": inexact,
                "cells_order": cells_order,
                "right_screw_angle": to_neighbor_angle_latex,
                # "tubes_to_3d": identify_tubes(Cs, v, P)
            }

            ## collect faces
            F_faces = []
            for face_equation in F.inequality_generator():
                face_vertices = cyclic_sort_vertices_2d([vert for vert in face_equation.incident()])
                vert_0, vert_1, vert_2 = [proj_sch(vector(v)) for v in face_vertices[:3]]
                face_normal = (vert_2 - vert_0).cross_product(vert_1 - vert_0)
                if face_normal.dot_product(vert_0) < 0:
                    face_vertices.reverse()
                F_faces.append([F.vertices().index(vert) for vert in face_vertices])
            F_vertices = [Q(p) for p in F.vertices_list()]
    
            ## collect vertices
            for pp in F_vertices:
                data["points"].append(tuple(proj_sch(pp.coefficient_tuple())))
            points_positions = {g:i for i, g in enumerate(data["points"])}

            ## convert to inexact quaternions for faster calculations
            Cs = [[FDR(Q_in(list(g.l)), g.r.to_quaternion()) for g in C] for C in Cs]
            F_vertices = [Q_in(list(pp)) for pp in F_vertices]


            ## construct tubes
            def point_index(gp):
                """
                return index of `gp` in `points_positions`.
                """
                for point, index in points_positions.items():
                    if math.sqrt(sum((point[ii]-gp[ii])**2 for ii in range(3))) <= ALMOST_ZERO:
                        return index
                return -1

            if verbose: print("TUBES: constructing... ", end=" ", flush=True)
            for ii in range(len(Cs)):
                visited_faces = set()
                C = Cs[ii]
                temp = data["tubes"].setdefault(str(ii), []);
                for g in C:
                    cell_indices = []
                    for pp in F_vertices:
                        gp = tuple(proj_sch(g.act(pp)))
                        try:
                            gp_index = points_positions[gp]
                        except KeyError:
                            gp_index = point_index(gp)
                        if gp_index == -1:
                            points_positions[gp] = len(data["points"])
                            gp_index = len(data["points"])
                            data["points"].append(gp)
                        cell_indices.append(gp_index)
                    cell_faces = []
                    for face in F_faces:
                        cell_face = [cell_indices[kk] for kk in face]
                        if frozenset(cell_face) in visited_faces: continue
                        visited_faces.add(frozenset(cell_face))
                        cell_faces.append(cell_face)
                    data["tubes"][str(ii)].append(cell_faces)
                if verbose: print(f"{ii}", end=" ", flush=True)
            if verbose: print("")
            ## turn point into tuples of type float
            ## this is needed so that json.dump works
            data["points"] = [tuple([float(_) for _ in p]) for p in data["points"]]

            ## write tube data
            output_directory = f"data/{gr_name_func(0)}/{k}"
            json_fname = f"{int(ncells)*int(ntubes)}cells_{ntubes}tubes.json"
            Path(output_directory).mkdir(parents=True, exist_ok=True)
            with open(f"{output_directory}/{json_fname}", "w") as f:
                #json.dump(data, f)
                text = json.dumps(data)
                f.write(compact_json(text))
            if verbose: print(f"TUBES: done in {time()-start_constructing}s")
            if verbose: print(f"## FINISH TUBE\n")

        ## write gallery data
        Path("data/cells/").mkdir(parents=True, exist_ok=True)
        json_fname = f"{gr_name_func(0)}_{k}.json"
        with open("data/cells/" + json_fname , "w") as f:
            json.dump(gallery, f)
        if verbose: print("## FINISH GROUP AND CENTER")

print(f"OVERALL RUNTIME: {time() - start_overall}")
