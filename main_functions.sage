#!/usr/bin/env python 
# coding: utf-8

## Contains definitions, classes and functions that are used in the
## other scripts.

from collections import Counter
from time import time
from pathlib import Path
import json
import sys

## quaternion_algebra
Q.<i,j,k> = QuaternionAlgebra(AA, -1, -1)
Q_in.<i_in,j_in,k_in> = QuaternionAlgebra(RDF, -1, -1)

## constants
ALMOST_ZERO = 1e-8

## useful unit quaternions [CS, Page 33]
sigma = (sqrt(5) - 1)/2
tau = (sqrt(5) + 1)/2
i_I = Q([0, 1/2, sigma/2, tau/2])
i_O = Q([0, 0, 1/sqrt(2), 1/sqrt(2)])
w = Q([-1/2, 1/2, 1/2, 1/2])
w_bar = w.conjugate()
e = lambda n: Q([cos(pi/n), sin(pi/n), 0, 0])

## other useful unit quaternions
sigma2 =  (-sqrt(5) - 1)/2
tau2 = (-sqrt(5) + 1)/2
i_I_dag = Q([0, 1/2, sigma2/2 , tau2/2]) # DuVal
i_I2 = Q([0, -sigma/2, -tau/2, 1/2]) # i'_I in [Cs, Table 4.2]

## useful functions
def nearby_rational_angle(angle):
    """
    Return a rational number `a` such that `a*pi` approximates `angle`.
    """
    return (angle/pi).n().nearby_rational(ALMOST_ZERO)

def quaternion_classify(el):
    """
    Classify the quaternion e by conjugacy type.
    """
    return nearby_rational_angle(arccos(el[0]))

class quat_2D:
    """
    Elements of the (binary) dihedral group 2D_2n.
    An instance `quat_2D(alpha, p)` represents the quaternion
    `(cos(alpha*pi) + i*sin(alpha*pi))*j^p`.
    Here 0<=alpha<2 and p is 0 or 1.
    """
    def __init__(self, alpha=0, p=0):
        self.alpha = quat_2D.mod2(alpha)
        self.p = p % 2

    def __eq__(self,other):
        return self.alpha == other.alpha and self.p == other.p

    def __hash__(self):
        return hash((self.alpha,self.p))

    def __mul__(self, other):
        if self.p == 0:
            return quat_2D(self.alpha+other.alpha, other.p)
        if other.p == 0:
            return quat_2D(self.alpha-other.alpha, 1)
        return quat_2D(1+self.alpha-other.alpha)

    def __pow__(self,n):
        if n == 0:
            return quat_2D(0)
        if n > 0:
            x=x0 = quat_2D(self.alpha,self.p)
        else:
            if self.p == 0:
                x=x0 = quat_2D(-self.alpha)
            else:
                x=x0 = quat_2D(1+self.alpha,-1)
        for _ in range(abs(n)-1):
            x *= x0
        return x

    def classify(self):
        """
        Classify by conjugacy type.
        """
        if self.p:
            return 1/2 # arccos(0)
        return min(self.alpha, 2-self.alpha)

    def conjugate(self):
        if self.p == 0:
            return quat_2D(2-self.alpha)
        else:
            return quat_2D(1+self.alpha, 1)

    def swap(self,n):
        """
        Conjugate (3d-rotation) self with the quaternion
        special_quaternions[n] (normalized).

        The following special quaternions:
        j+k (n=4), j-k (n=5), 1+i (n=6), i(n=9), j(n=10), k(n=11)
        normalize the group Dn_group.
        The other special quaternions only normalize the group
        Q8 = {+-1, +-i, +-j, +-k}.
        """
        # identity
        if n<0: return self
        
        special_quaternion = special_quaternions[n]

        # ones which preserve Dn_group
        if special_quaternion == j+k:
            if self.p: return quat_2D(1/2-self.alpha, 1)
            else: return quat_2D(-self.alpha, 0)
        if special_quaternion == j-k:
            if self.p: return quat_2D(3/2-self.alpha, 1)
            else: quat_2D(-self.alpha, 0)
        if special_quaternion == 1+i:
            if self.p: return quat_2D(self.alpha-1/2, 1)
            else: return quat_2D(self.alpha, 0)
        if special_quaternion == i:
            if self.p: return quat_2D(1+self.alpha, 1)
            else: return quat_2D(self.alpha, 0)
        if special_quaternion == j:
            if self.p: return quat_2D(-self.alpha, 1)
            else: return quat_2D(-self.alpha, 0)
        if special_quaternion == k:
            if self.p: return quat_2D(1-self.alpha, 1)
            else: return quat_2D(-self.alpha, 0)
        
        # one which preserve Q8
        if self.alpha not in (0, 1/2, 1, 3/2):
            raise ValueError("alpha not a multiple of 1/2")
        return self.swap_dict[n][self]
        
    def mod2(alpha):
        while alpha < 0:
            alpha += 2
        while alpha >= 2:
            alpha -= 2
        return alpha

    def to_quaternion(self):
        return (cos(self.alpha*pi).n()+i_in*sin(self.alpha*pi).n()) * j_in^self.p

    def to_exact_quaternion(self):
        return (AA(cos(self.alpha*pi)) + i*AA(sin(self.alpha*pi))) * j^self.p

    def __str__(self):
        if self.p:
            if self.alpha == 0:
                return "j"
            elif self.alpha == 1:
                return "-j"
            elif self.alpha == 1/2:
                return "k"
            elif self.alpha == 3/2:
                return "-k"
            return f"e^{self.alpha}*j"
        else:
            if self.alpha==0:
                return "1"
            elif self.alpha==1:
                return "-1"
            elif self.alpha==1/2:
                return "i"
            elif self.alpha==3/2:
                return "-i"
            return f"e^{self.alpha}"

    def __repr__(self):
        return f"quat_2D({self.alpha},{self.p})"

## usually used quat_2D elements
quat_i = quat_2D(1/2,0)
quat_j = quat_2D(0,1)
quat_k = quat_i*quat_j
quat_1 = quat_2D(0)
quat_minus1 = quat_2D(1)

## fill quat_2D.swap_dict
special_quaternions = [i+j, i-j, i+k, i-k, j+k, j-k,
                       1+i, 1-i, 1+j, 1-j, 1+k, 1-k,
                       i, j, k,
                       1+i+j+k, 1-i+j+k, 1+i-j+k, 1+i+j-k,
                       1-i-j+k, 1-i+j-k, 1+i-j-k, 1-i-j-k
                    ]
def fill_swap_dict(verbose=False):
    units = list(a*b for a in (quat_1, quat_minus1) for b in (quat_1, quat_i, quat_j, quat_k))
    quat_2D.swap_dict = []
    for n,q in enumerate(special_quaternions):
        if verbose:
            for e in (i,j,k):
                print (f"type {n}. conjugation with ({q}) maps",e,"->",q^-1*e*q)
        res = dict()
        q = Q_in(list(q))
        for e in units:
            ee = q^-1*e.to_quaternion()*q
            for f in units:
                if norm(vector(f.to_quaternion()-ee))<1e-2:
                    res[e]=f
                    break
            else:
                raise "something is wrong"
        quat_2D.swap_dict.append(res)
    return
fill_swap_dict()

def tests_on_swap_dict():
    units = list(a*b for a in (quat_1, quat_minus1) for b in (quat_1, quat_i, quat_j, quat_k))
    for n,q in enumerate(special_quaternions):
        for e in units:
            ee = q^-1*e.to_exact_quaternion()*q
            assert ee == e.swap(n).to_exact_quaternion(), n
    return
# tests_on_swap_dict()

class FDR: 
    """
    A 4-dimensional orthogonal transformation [l, r] or *[l, r], where
    [l, r]x: x -> l.conjugate() * x * l
    is orientation preserving and,
    *[l, r]x: x -> l.conjugate() * x.conjugate() * l
    is orientation reversing.

    INPUT:
    - `l` -- quaternion or ``quat_2d``: A unit quaternion.
    - `r` -- quaternion or ``quat_2d``: A unit quaternion.
    - `star` -- boolean (default: ``False``): If ``False`` construct [l, r].
       If `True` construct *[l, r].
    """
    def __init__(self, l, r, star=False):
        self.l = l
        self.r = r
        self.star = star

    def act(self, x):
        """
        Return that action of self on x, where x is a quaternion.
        """
        if self.star:
            x = x.conjugate()

        l_bar = self.l.conjugate()
        if type(l_bar) is quat_2D:
            l_bar = l_bar.to_exact_quaternion()
        r = self.r
        if type(r) is quat_2D:
            r = r.to_exact_quaternion()
        return l_bar*x*r

    def __mul__(self,other):
        if other.star:
            return FDR(self.r*other.l, self.l*other.r, not self.star)
        else:
            return FDR(self.l*other.l, self.r*other.r, self.star)

    def inverse(self):
        if self.star:
            return FDR(self.r.conjugate(), self.l.conjugate(), True)
        else:
            return FDR(self.l.conjugate(), self.r.conjugate())

    def classify(self):
        """
        Classify self by its conjugacy type.

        For [l, r] return (a, b) with 0 <= a, b <= 1 where
        [l, 1] is a left rotation by +/- a*pi and
        [1, r] is a right rotation by +- b*pi.
        The pair [-l,-r], which represents the same rotation,
        will yield the pair (1-a,1-b). 
        (Since the groups that are considered
        are closed under taking negatives,
        both equivalent possibilities will be present.)
        Hence, we normalize by requiring that a<b or a=b<=1/2.
        
        For *[l, r] return (7, b) with 0 <= b <= 1
        representing a rotation-reflection by +/- (1-b)*pi.
        """
        if self.star:
            try:
                lr = self.l*self.r
            # except KeyError:
            except (TypeError, AttributeError) as e:
                if type(self.l) is quat_2D:
                    lr = self.l.to_exact_quaternion()*self.r
                else:
                    lr = self.l * self.r.to_exact_quaternion()
            try:
                c = lr.classify()
            except AttributeError:
                c = quaternion_classify(lr)
            return (7,c)
            
        try:
            cl = self.l.classify()
        except AttributeError:
            cl = quaternion_classify(self.l)
        try:
            cr = self.r.classify()
        except AttributeError:
            cr = quaternion_classify(self.r)
        if cl > cr or cl == cr > 1/2:
            return (1-cl,1-cr) # normalize
        return (cl,cr)

    def __str__(self):
        brak = f"[{str(self.l)}, {str(self.r)}]"
        if self.star:
            return "*"+brak
        return brak

    def __repr__(self):
        s = ""
        if self.star:
            s = ",1"
        return f"FDR({self.l!r}, {self.r!r}{s})"

    def __eq__(self,other):
        return self.l == other.l and self.r == other.r and self.star == other.star

    def __hash__(self):
        return hash((self.l, self.r, self.star))

    def to_exact_quaternions(self):
        """
        return self with quaternion left and right components.
        """
        newl = self.l.to_exact_quaternion() if isinstance(self.l, quat_2D) else self.l
        newr = self.r.to_exact_quaternion() if isinstance(self.r, quat_2D) else self.r
        return FDR(newl, newr, self.star)

    def to_matrix(self):
        """
        return matrix represnetation of self.
        """
        M = [(self.act(a)).coefficient_tuple() for a in Q.basis()]
        return matrix(Q.base_ring(), M).transpose()

## mixed one
one_HQ = FDR(Q(1),quat_1)
one_QH = FDR(quat_1,Q(1))

def element_chars(el):
    """
    If `el` is instance of `quat_2D`, return `(quat_1, quat_minus1)`.
    Else, return `(Q(1), Q(-1))`.
    """
    if isinstance(el, quat_2D):
        return quat_1, quat_minus1#, Dn_group # or Cn_group
    return Q(1), Q(-1)#, Q_group

def mk_group(gens, one=None, limit=30000):
    """
    Return the group generated by elements in `gens`.
    
    INPUT:
    - `gens` -- list: Generators.
    - `one` -- any: The one of the group.
    - `limit` -- int: Run time bound on the order of the group.
    """
    if one is None:
        if isinstance(gens[0], FDR):
            one = FDR(element_chars(gens[0].l)[0], element_chars(gens[0].r)[0])
        else:
            one = Q(1)
    # maintain group as a list because group 
    # elements are generated in a unique order
    G = [one] 
    G_set = {G[0]}
    new_in_G = list(G)
    while new_in_G:
        to_come = []
        for g in new_in_G:
            for h in gens:
                new_g = g*h
                if new_g not in G_set:
                    G.append(new_g)
                    G_set.add(new_g)
                    to_come.append(new_g)
                    if limit and len(G) > limit:
                        raise ValueError("group size over the limit")
        new_in_G = to_come
    return G

## Hurley
Hurley_table = { # [Brown et al., Page 55]
    (-4, 6, 1): ("I'",'2222'),
    (-3, 4, 1): ("Z'",'322'),
    (-2, 0, -1): ("T'",'2221'),
    (-2, 2, 1): ("R'",'422'),
    (-2, 3, 1): ("S'",'33'),
    (-1, 0, -1): ("N'",'321'),
    (-1, 0, 1): ("K'",'622'),
    (-1, 1, 1): ("L'",'5'),
    (-1, 2, 1): ("M'",'43'),
    (0, -2, 1): ("E",'2211'),
    (0, -1, 1): ("C",'T=12'),
    (0, 0, -1): ("F",'421'),
    (0, 0, 1): ("A",'8'),
    (0, 1, 1): ("B",'63'),
    (0, 2, 1): ("D",'44'),
    (1, 0, -1): ("N",'621'),
    (1, 0, 1): ("K",'311'),
    (1, 1, 1): ("L",'X=10'),
    (1, 2, 1): ("M",'64'),
    (2, 0, -1): ("T",'2111'),
    (2, 2, 1): ("R",'411'),
    (2, 3, 1): ("S",'66'),
    (3, 4, 1): ("Z",'611'),
    (4, 6, 1): ("I",'1111'),
}

Hurley_conversion = {
    (7, 1/3): (1, 0, -1),
    (11/12, 7/12): (1, 2, 1),
    (1/6, 1/2): (0, 1, 1),
    (1/2, 5/6): (0, 1, 1),
    (7, 2/3): (-1, 0, -1),
    (1, 2/3): (2, 3, 1),
    (1/12, 5/12): (1, 2, 1),
    (1/2, 1/2): (0, -2, 1),
    (1/4, 3/4): (-2, 2, 1),
    (1/2, 1/3): (0, -1, 1),
    (2/3, 0): (-2, 3, 1),
    (1/3, 0): (2, 3, 1),
    (5/12, 1/12): (1, 2, 1),
    (2/3, 1): (2, 3, 1),
    (1/3, 1): (-2, 3, 1),
    (1/6, 1/6): (3, 4, 1),
    (2/3, 1/3): (-1, 0, 1),
    (3/4, 3/4): (2, 2, 1),
    (1/3, 1/3): (1, 0, 1),
    (1/2, 3/4): (0, 0, 1),
    (1, 0): (-4, 6, 1),
    (7/12, 1/12): (-1, 2, 1),
    (5/6, 5/6): (3, 4, 1),
    (1/5, 3/5): (-1, 1, 1),
    (0, 0): (4, 6, 1),
    (2/3, 2/3): (1, 0, 1),
    (1/3, 2/3): (-1, 0, 1),
    (0, 1): (-4, 6, 1),
    (5/6, 1/2): (0, 1, 1),
    (1/4, 1/2): (0, 0, 1),
    (11/12, 5/12): (-1, 2, 1),
    (7, 1/2): (0, 0, -1),
    (4/5, 2/5): (-1, 1, 1),
    (1, 1/3): (-2, 3, 1),
    (0, 1/3): (2, 3, 1),
    (1/4, 1/4): (2, 2, 1),
    (4/5, 3/5): (1, 1, 1),
    (0, 2/3): (-2, 3, 1),
    (1/2, 1): (0, 2, 1),
    (3/4, 1/2): (0, 0, 1),
    (5/6, 1/6): (-3, 4, 1),
    (3/4, 1/4): (-2, 2, 1),
    (1/12, 7/12): (-1, 2, 1),
    (1/2, 1/4): (0, 0, 1),
    (2/5, 1/5): (1, 1, 1),
    (3/5, 1/5): (-1, 1, 1),
    (1/2, 2/3): (0, -1, 1),
    (2/3, 1/2): (0, -1, 1),
    (1/5, 2/5): (1, 1, 1),
    (1/3, 1/2): (0, -1, 1),
    (5/12, 11/12): (-1, 2, 1),
    (1/2, 1/6): (0, 1, 1),
    (7, 0): (2, 0, -1),
    (7, 1): (-2, 0, -1),
    (2/5, 4/5): (-1, 1, 1),
    (1, 1): (4, 6, 1),
    (3/5, 4/5): (1, 1, 1),
    (1, 1/2): (0, 2, 1),
    (7/12, 11/12): (1, 2, 1),
    (0, 1/2): (0, 2, 1),
    (1/2, 0): (0, 2, 1),
    (1/6, 5/6): (-3, 4, 1)
}

def Hurley_pattern(pattern):
    """
    Convert "our" classification code (pairs (a,b) with multiplicities)
    into Hurley-pattern; divide multiplicities by 2
    """
    try:
        H_cla = {}
        for code in pattern.split():
            # decode into class,num
            class_code, num_code = code.split(":")
            num = Integer(num_code)
            if class_code[0]=="*":
                clas = 7, Rational(class_code[1:])
            else:
                if "/" in class_code:
                    a, denominator_code = class_code.split("/")
                    denominator = Integer(denominator_code)
                    class_code = a
                else:
                    denominator = 1
                a, b = class_code.split("|")
                clas = Integer(a)/denominator, Integer(b)/denominator
            h3 = Hurley_conversion[clas]
            hcode, _ = Hurley_table[h3]
            H_cla[hcode] = H_cla.get(hcode,0) + num

        H_cl = []
        for k,num in sorted(H_cla.items()):
            assert num % 2 == 0
            H_cl.append(f"{num//2}*{k}")
        return " ".join(H_cl)

    except KeyError:
        return # definitely not a crystallographic group.

# print("Hermann-symbol / our class. / Hurley symbol / Hurley triplet")
# for x in sorted((he,cl,hu,hurley3) for cl,hurley3 in Hurley_conversion.items()
#     for hu,he in (Hurley_table[hurley3],)):
#         print("{0:6} ({1[0]!s:5},{1[1]!s:5})   {2:4} ({3[0]:-2},{3[1]:3},{3[2]:3})".format(*x))


class Q_group:
    """
    Quaternion group.

    INPUT:
    - `gens` -- list of quaternions or `quat_2D`'s: Generators of the group.
    - `limit` -- int: Run time bound on the order of the group.
    """
    def __init__(self, gens, limit=None):
        self.gens = gens
        self.one, self.minus_1 = element_chars(gens[0])
        self.q = mk_group(self.gens, self.one, limit)
        self.order = len(self.q)

    def ensure_number(self):
        """
        Construct `self.p_conjugate`, the permutation that maps
        every element to its conjugate (=inverse).
        """
        if hasattr(self,"number"):
            return
        self.number = {e:i for i,e in enumerate(self.q)} # start from 0
        self.p_conjugate = Permutation([self.number[e.conjugate()]+1 for e in self.q])

    def ensure_pgroup(self):
        """
        construct `self.p_group`, the permutation group of `self`.
        """
        if hasattr(self,"pgroup"):
            return
        self.ensure_number()
        self.pgens = [self.mk_permutation(e) for e in self.gens]
        
        pgens_cycles = [Permutation([x+1 for x in perm]).cycle_string() for perm in self.pgens]
        self.pgroup = PermutationGroup(pgens_cycles)

    def ensure_name_from_perm(self):
        if hasattr(self,"name_from_perm"):
            return
        self.name_from_perm = {tuple(self.mk_permutation(e)):e for e in self.q}
        
    def mk_permutation(self,e):
        """
        right multiplication by e, represented as a permutation.
        """
        return ([self.number[g*e] for g in self.q])
    
class Dn_group(Q_group):
    """
    Quaternion group 2D_n (`n` even) of order 2n.
    """
    def __init__(self, n):
        assert n % 2 == 0       
        gens = [quat_2D(2/n), quat_2D(0,1)]
        Q_group.__init__(self, gens)

class Cn_group(Q_group):
    """
    Quaternion group 2C_n of order 2n.
    """
    def __init__(self, n):
        gens = [quat_2D(1/n)]  
        Q_group.__init__(self, gens)


def compress_pair(a, b):
    """
    Compressed notation for the rational pair `(a, b)`. This is,
    `a*n|b*n/n where n is the common denominator of `a` and `b`.
    """
    if a == 7:
        return f"*{b!s}"
    an,ad = a.numerator(),a.denominator()
    bn,bd = b.numerator(),b.denominator()
    n = lcm(ad,bd)
    res = f"{a*n}|{b*n}"
    if n>1:
        res += f"/{n}"
    return res

def classify_group(gr):
    """
    Classify a group according to the multiplicities of its
    elements in the different geometric conjugacy classes.
    `gr` is a list of group elements.
    """
    multi = Counter()
    for el in gr:
        c = el.classify()
        multi[c] += 1
    # turn classification into a compact string
    return " ".join(f"{compress_pair(a,b)}:{c}" for (a,b),c in sorted(multi.items()))

def group_string_rep(gr):
    """
    A way to represent the group as a string.
    """
    if type(gr) != list:
        gr = gr.group
    return " ".join(sorted(map(str,gr)))

class AxB_group:
    """
    AxB group of a 4-dimensional point group.

    In this group, both [l, r] and [-l, -r] are presented as distinct
    elements. Hence, this group contains twice the number of elements of
    the corresponding point group.

    INPUT:
    - `gens` -- list of FDR's: Generators of the group.
    - `name` -- str: Group name.
    - `limit` -- int: Run time bound on the order of the group.
    - `group` -- list of FDR's: In case group is already generated.
      In this case, gens is also required and should has the
      same data type as group.
    """
    def __init__(self, gens, name=None, limit=None, group=None):
        self.gens = gens
        self.name = name
        self.limit = limit
        self.achiral = any(e.star for e in gens)
        l_one, l_minus1 = element_chars(gens[0].l)
        r_one, r_minus1 = element_chars(gens[0].r)
        self.one = FDR(l_one,r_one)
        self.minus_1_minus_1 = FDR(l_minus1, r_minus1)
        self.group = group if group else mk_group(gens, one=self.one, limit=limit)
        self.order = len(self.group)

    def show(self):
        return list(str(e) for e in self.group)

    def __repr__(self):
        name = "a point group"
        if self.name: name = self.name
        return f"AxB group of {name} of order {self.order//2}"
 
    def __eq__(self, other):
        return set(self.group) == set(other.group)
                 
    def classify_group(self):
        """
        Classify the group according to the multiplicities of
        its elements in the conjugacy classes.
        """
        return classify_group(self.group)

    def Hurley_pattern(self):
        """
        Classify the group according the conjugacy classes of
        Hurley for crystallographic groups.
        """
        return Hurley_pattern(classify_group(self.group))

    def ensure_groups(self):
        """
        Generate the left and the right quaternion groups Q1 and Q2.
        """
        if hasattr(self,"Q1") and hasattr(self.Q1,"number"):
            return        
        self.Lgens = [e.l for e in self.gens]
        self.Rgens = [e.r for e in self.gens]
        if self.achiral:
            self.Lgens += self.Rgens
            self.Rgens = self.Lgens
            self.Q1 = self.Q2 = Q_group(self.Lgens)
            self.Q1.ensure_number()
        else:
            self.Q1 = Q_group(self.Lgens)
            self.Q2 = Q_group(self.Rgens)
            self.Q1.ensure_number()
            self.Q2.ensure_number()
        self.nL = self.Q1.order
        self.nR = self.Q2.order

    def ensure_pgroup(self):
        """
        If self is chiral, construct the permutation group of the AxB group.
        If self is achiral, construct the permutation group of the
        corresponding point group. Hence, its order is half the order
        of the AxB group.
        """
        if hasattr(self,"pgroup"):
            return
        self.ensure_groups()
        if not self.achiral:
            pgens = []
            for g in self.gens:
                pl = self.Q1.mk_permutation(g.l)
                pr = self.Q2.mk_permutation(g.r)
                pair = self.pair_permutation(pl,pr)
                pgens.append(pair)
            self.pgroup = PermutationGroup(pgens)
            pl = self.Q1.mk_permutation(self.Q1.minus_1)
            pr = self.Q2.mk_permutation(self.Q2.minus_1)
            self.Pminus_1_minus_1 = self.pair_permutation(pl,pr)
            assert self.pgroup.order() == self.order, (self.name,"WR-CHIRAL")
            return

        else: # achiral
            # represent as a permutation group acting on some elements.
            # Add i,j,k to the elements, to make sure the action is faithful.
            # But [-1,-1] is the identity permutation. So we only get HALF
            # of the elements.
            # This would also be possible for the chiral groups, but
            # then it would not be easy to translate the permutations
            # into [l,r] pairs.

            def to_permutation(g):
                """permutation as sequence of numbers 1,2,...."""
                if g.star:
                    return (elementnumber[
                        g.l.conjugate() * e1.conjugate() * g.r ]+1
                         for e1 in elements)
                else:
                    return (elementnumber[
                        g.l.conjugate() * e1 * g.r ]+1 for e1 in elements)
            
            if type(self.group[0].l) is quat_2D:
                extension_elements = [quat_i,quat_j,quat_k]
            else:
                extension_elements = [i,j,k]

            elements = mk_group(self.Q1.gens + extension_elements, self.Q1.one, self.limit)
            elementnumber = {e:i for i,e in enumerate(elements)}
            # start from 0
                
            self.to_AxB_group = {tuple(to_permutation(g)): g for g in self.group}
            pgens = [Permutation(to_permutation(g)) for g in self.gens]
            self.pgroup = PermutationGroup(pgens)
            assert 2*self.pgroup.order() == self.order, (self.name,"WR-ACHIRAL")
        
    def pair_permutation(self, perm1, perm2): 
        """
        Representation of [l,r] as a permutation, where
        `l=perm1` and `r=perm2` are represented as permutations
        (lists starting at 0).
        """
        assert not self.achiral
        p = Permutation([x+1 for x in perm1]+[x+self.nL+1 for x in perm2])
        return p.cycle_string()
    
    def is_isomorphic(self,other):
        self.ensure_pgroup()
        other.ensure_pgroup()
        return self.pgroup.is_isomorphic(other.pgroup)

    def conjugate_group(self, q1, q2=None, scale1=None, scale2=None):
        """
        Conjugate the group elements.
        See self.conjugate_element documentation.
        """
        return [self.conjugate_element(f, q1, q2, scale1, scale2) for f in self.group]
        
    @staticmethod
    def conjugate_element(f, q1, q2=None, scale1=None, scale2=None):
        """
        Conjugate `f` by `FDR(q1, q2)`.
        Here `q1` (and `q2`) can be either:
        * quaternion: if it is not normalized, then its norm squared,
            scale1, should be provided.
        * a 2-tuple `("S", n)`: meaning `f.l.swap(n)`.
        * a 3-tuple `("S", n, qq1)`: meaning `qq1.conjugate()*f.l.swap(n)*qq1`.
        """
        # set up left and right conjugation maps.
        if type(q1) is tuple and q1[0]=="S":
            n1 = q1[1]
            if f.star:
                if q1==q2 and len(q1)==2:
                    return FDR(f.l.swap(n1),f.r.swap(n1),True)
                q2,scale2 = q1,scale1 #
                raise NotImplementedError
            if len(q1)==2:
                conjugate_left = lambda e: e.swap(n1)
            else:
                qq1 = q1[2]
                conjugate_left = lambda e: qq1.conjugate()*e.swap(n1)*qq1
                # could try swap before or after
        elif q1 is None:
            conjugate_left = lambda e: e
        elif scale1:
            conjugate_left = lambda e: q1.conjugate()*e*q1*scale1
        else:
            conjugate_left = lambda e: q1.conjugate()*e*q1

        if type(q2) is tuple and q2[0]=="S":
            if f.star:
                raise NotImplementedError
            n2 = q2[1]
            if len(q2)==2:
                conjugate_right = lambda e: e.swap(n2)
            else:
                qq2 = q2[2]
                conjugate_right = lambda e: qq2.conjugate()*e.swap(n2)*qq2
        elif q2 is None:
            conjugate_right = lambda e: e
        elif scale2:
            conjugate_right = lambda e: q2.conjugate()*e*q2*scale2
        else:
            conjugate_right = lambda e: q2.conjugate()*e*q2
            
        if f.star:
            if type(f.l) is quat_2D:
                if scale1 not in [None, quat_1]:
                    raise ValueError
                if scale2 not in [None, quat_1]:
                    raise ValueError
                scale1s = quat_1
                scale2s = quat_1
            else:
                scale1s = AA(sqrt(scale1)) if scale1 else 1
                scale2s = AA(sqrt(scale2)) if scale2 else 1
            
            if q2 is None:
                return FDR(f.l*q1*scale1s, q1.conjugate()*f.r*scale1s, True)
            else:
                return FDR(q2.conjugate()*f.l*q1*scale1s*scale2s, q1.conjugate()*f.r*q2*scale1s*scale2s, True)
        else:
            return FDR(conjugate_left(f.l), conjugate_right(f.r), False)


    def find_conjugation(self, other):
        """
        Try to find conjugation between the two AxB groups.
        If the groups are not conjugate, return `False`.
        If the groups are conjugate, try to find the conjugation.
        If the conjugation if found, return `True` and set the values
        for `self.conj_name` and `self.conj_latex`.
        If the conjugation is not found, raise `NotImplementedError`.

        IMPORTANT:
        The procedure does not guarantee to find the conjugation. It
        only tests a set of special conjugations.
        If the components of the elements of self are `quat_2D`, then
        the procedure is faster. However, it tries a smaller set of 
        special quaternions.
        If the components of the elements of self are quaternions,
        the procedure tries a larger set of special conjugations.
        Thus, it is recommended to convert the elements to be so
        by using `self.to_exact_quaternion_group()`. However, this is
        slower.
        If the components of the elements of self are mixed of quaternions
        and `quat_2D`'s, then the procedure starts by converting both
        components to be quaternions by calling 
        `self.to_exact_quaternion_group()`.
        """
        def check_equality(gr):
            """
            Check if gr and other are the same. 
            """
            if group_string_rep(gr) == other_string:
                return True
            return False
       
        def check_equality_elementwise(gr):
            """
            Check if gr and other (H) are the same.
            (This is used as a fallback if check_equality fails)
            (Only used if components of elements of gr are quaternions)
            """
            if set(gr) == H_set:
                return True
            return False

        def swap_to_latex(n, qq=Q(1)):
            """
            Convert the given conjugation to latex
            """
            if isinstance(qq, quat_2D):
                qq = qq.to_exact_quaternion()
            if n < 0:
                q = qq
            else:
                q = special_quaternions[n]*qq
            if q.reduced_norm() == 2:
                s = "\\frac1\\sqrt2"
            elif q.reduced_norm() == 4:
                s = "\\frac12"
            elif q.reduced_norm() == 1:
                s = "1"
            else:
                raise 
            return f"{s}({q})"

        ## check if they are conjugate
        if self.classify_group() != other.classify_group():
            return False

        ## construct group string of other
        other_string = group_string_rep(other.group)
        assert other_string == " ".join(sorted(other.show()))

        ## conjugations for quat_2D groups
        only_quat_2D = type(self.group[0].l) is quat_2D and type(self.group[0].r) is quat_2D
        if only_quat_2D:

            ## conjugation by elements of Q8
            units = list(a*b for a in (quat_1, quat_minus1) for b in (quat_1, quat_i, quat_j, quat_k))
            for unit1 in units:
                for unit2 in units:
                    G_conj = self.conjugate_group(unit1, unit2)
                    if check_equality(G_conj):
                        self.conj_name = f"[{unit1}, {unit2}]"
                        self.conj_latex = f"[{unit1}, {unit2}]"
                        return True

            ## S conjugations from left
            for n in range(-1, len(special_quaternions)):
                try:
                    G_conj = self.conjugate_group(("S", n))
                    if check_equality(G_conj):
                        self.conj_name = f"[S_{{{n}}}, 1]"
                        self.conj_latex = f"[{swap_to_latex(n)}, 1]"
                        return True
                    # for unit in [quat_i, quat_j, quat_k]:
                    #     G_conj_temp = [self.conjugate_element(e, unit) for e in G_conj]
                    #     if check_equality(G_conj_temp):
                    #         self.conj_name = f"[S_{{{n}}}, 1][{unit.to_exact_quaternion()}, 1]"
                    #         self.conj_latex = f"[{swap_to_latex(n, unit)}, 1]"
                    #         return True
                except (ValueError, NotImplementedError):
                    continue

            ## S conjugations from right
            for n in range(0, len(special_quaternions)):
                try:
                    G_conj = self.conjugate_group(("S", -1), ("S", n))
                    if check_equality(G_conj):
                        self.conj_name = f"[1, S_{{{n}}}]"
                        self.conj_latex = f"[1, {swap_to_latex(n)}]"
                        return True
                    # for unit in [quat_i, quat_j, quat_k]:
                    #     G_conj_temp = [self.conjugate_element(e, unit) for e in G_conj]
                    #     if check_equality(G_conj_temp):
                    #         self.conj_name = f"[1, S_{{{n}}}][{unit.to_exact_quaternion()}, 1]"
                    #         self.conj_latex = f"[1, {swap_to_latex(n, unit)}]"
                    #         return True
                except (ValueError, NotImplementedError):
                    continue

            ## S conjugations from both sides
            for n1 in range(0, len(special_quaternions)):
                for n2 in range(0, len(special_quaternions)):
                    try:
                        G_conj = self.conjugate_group(("S", n1), ("S", n2))
                        if check_equality(G_conj):
                            self.conj_name = f"[S_{{{n1}}}, S_{{{n2}}}]"
                            self.conj_latex = f"[{swap_to_latex(n1)}, {swap_to_latex(n2)}]"
                            return True
                        # for unit1 in [quat_i, quat_j, quat_k]:
                        #     for unit2 in [quat_i, quat_j, quat_k]:
                        #         G_conj_temp = [self.conjugate_element(e, unit1, unit2) for e in G_conj]
                        #         if check_equality(G_conj_temp):
                        #             self.conj_name = f"[S_{{{n1}}}, S_{{{n2}}}][{unit1.to_exact_quaternion()}, {unit2.to_exact_quaternion()}]"
                        #             self.conj_latex = f"[{swap_to_latex(n1, unit1)}, {swap_to_latex(n2, unit2)}]"
                        #             return True
                    except (ValueError, NotImplementedError):
                        continue

            # (S, alpha) conjugations from both sides       
            # look for elements of the form [(alpha,1),r] or [r,(alpha,1)] without star
            alpha_left2  = set(e.l.alpha for e in other.group if not e.star and e.l.p)
            alpha_right2 = set(e.r.alpha for e in other.group if not e.star and e.r.p)
            for n1 in range(-1, len(special_quaternions)):
                for n2 in range(-1, len(special_quaternions)):
                    try:
                        G_conj = self.conjugate_group(("S", n1), ("S", n2))
                        alpha_left  = set(e.l.alpha for e in G_conj if not e.star and e.l.p)
                        alpha_right = set(e.r.alpha for e in G_conj if not e.star and e.r.p)
                        diffs_left  = set((alpha-alpha2)/2 for alpha in alpha_left for alpha2 in alpha_left2)
                        diffs_right = set((alpha-alpha2)/2 for alpha in alpha_right for alpha2 in alpha_right2)
                        if not diffs_left: diffs_left = [0]
                        if not diffs_right: diffs_right = [0]
                        diffs_left = self.sort_by_abs_value(diffs_left)
                        diffs_right = self.sort_by_abs_value(diffs_right)

                        for dl in diffs_left:
                            for dr in diffs_right:
                                try:
                                   G_conj_temp = [self.conjugate_element(e, quat_2D(dl), quat_2D(dr)) for e in G_conj]
                                   if check_equality(G_conj_temp):
                                       self.conj_name = f"[S_{{{n1}}}, S_{{{n2}}}][{quat_2D(dl)}, {quat_2D(dr)}]"
                                       self.conj_latex = f"[{swap_to_latex(n1, quat_2D(dl))}, {swap_to_latex(n2, quat_2D(dr))}]"
                                       return True
                                except (ValueError, NotImplementedError):
                                    continue
                    except (ValueError, NotImplementedError):
                        continue

            ## (S,q) conjugations from both sides 
            ## seems to be not needed
            if False:
                self.ensure_groups()
                for n1 in range(-1, len(special_quaternions)):
                    for q1 in self.Q1.q:
                        for n2 in range(0, len(special_quaternions)):
                            for q2 in self.Q2.q:
                                try:
                                    G_conj = self.conjugate_group(("S", n1, q1), ("S", n2, q2))
                                    if check_equality(G_conj):
                                        self.conj_name = f"[S_{{{n1}}}, S_{{{n2}}}][{q1}, {q2}]"
                                        self.conj_latex = f"[{swap_to_latex(n1, q1)}, {swap_to_latex(n2, q2)}]"
                                        return True
                                except (ValueError, NotImplementedError):
                                    continue
        
        ## conjugations for quaternion groups
        G = self.to_exact_quaternion_group()
        H = other.to_exact_quaternion_group()
        H_set = set(H.group)

        ## add square norm to special quaternions
        special_quaternions_snorms = (
                [(x, 1/2) for x in (i+j, i-j, i+k, i-k, j+k, j-k)] + 
                [(x, 1/2) for x in (1+i, 1-i, 1+j, 1-j, 1+k, 1-k)] +
                [(x, None) for x in (i, j, k)] +
                [(x, 1/4) for x in (1+i+j+k, 1-i+j+k, 1+i-j+k, 1+i+j-k)] +
                [(x, 1/4) for x in (1-i-j+k, 1-i+j-k, 1+i-j-k, 1-i-j-k)] +
                [(x, 1) for x in (e(8), e(-8))]
                )

        ## special conjugations from left
        for q,s in special_quaternions_snorms:
            G_conj = G.conjugate_group(q, None, s)
            if check_equality(G_conj):
                self.conj_name = f"[{q}, 1]"
                self.conj_latex = f"[{swap_to_latex(-1, q)}, 1]"
                return True
            if check_equality_elementwise(G_conj):
                self.conj_name = f"[{q}, 1]"
                self.conj_latex = f"[{swap_to_latex(-1, q)}, 1]"
                return True

        ## special conjugations from right
        for q,s in special_quaternions_snorms:
            G_conj = G.conjugate_group(Q(1), q, None, s)
            if check_equality(G_conj):
                self.conj_name = f"[1, {q}]"
                self.conj_latex = f"[1, {swap_to_latex(-1, q)}]"
                return True
            if check_equality_elementwise(G_conj):
                self.conj_name = f"[1, {q}]"
                self.conj_latex = f"[1, {swap_to_latex(-1, q)}]"
                return True

        ## special conjugations from both sides
        for q1,s1 in special_quaternions_snorms:
            for q2,s2 in special_quaternions_snorms:
                G_conj = G.conjugate_group(q1, q2, s1, s2)
                if check_equality(G_conj):
                    self.conj_name = f"[{q1}, {q2}]"
                    self.conj_latex = f"[{swap_to_latex(-1, q1)}, {swap_to_latex(-1, q2)}]"
                    return True
                if check_equality_elementwise(G_conj):
                    self.conj_name = f"[{q1}, {q2}]"
                    self.conj_latex = f"[{swap_to_latex(-1, q1)}, {swap_to_latex(-1, q2)}]"
                    return True
         
        ## group conjugations from both sides
        G.ensure_groups()
        for q1 in G.Q1.q:
            for q2 in G.Q2.q:
                G_conj = G.conjugate_group(q1, q2)
                if check_equality(G_conj):
                    self.conj_name = f"[{q1}, {q2}]"
                    self.conj_latex = f"[{swap_to_latex(-1, q1)}, {swap_to_latex(-1, q2)}]"
                    return True
                if check_equality_elementwise(G_conj):
                    self.conj_name = f"[{q1}, {q2}]"
                    self.conj_latex = f"[{swap_to_latex(-1, q1)}, {swap_to_latex(-1, q2)}]"
                    return True
        
        ## did not find the conjugation
        raise NotImplementedError
        return
    
    @staticmethod
    def sort_by_abs_value(s):
        return list(-v for _,v in sorted((abs(v0),-v0) for v0 in s))

    def to_exact_quaternion_group(self):
        """
        Return `self` as AxB group with both components of the elements
        being quaternions.
        """
        if isinstance(self.group[0].l, type(Q([0,1,0,0]))) and isinstance(self.group[0].r, type(Q([0,1,0,0]))):
            return self
        convert = lambda x: x.to_exact_quaternion() if isinstance(x, quat_2D) else x
        gens = [FDR(convert(e.l), convert(e.r), e.star) for e in self.gens]
        group = [FDR(convert(e.l), convert(e.r), e.star) for e in self.group]
        return AxB_group(gens, name=self.name, limit=self.limit, group=group)

    def split(self, e):
        """
        Split a permutation into two lists of tuples
        """
        perm=e.cycle_tuples()
        #print (perm)
        pleft = [c for c in perm if all(x <=self.nL for x in c)]
        pright = [c for c in perm if all(x >self.nL for x in c)]
        assert len(pleft) + len(pright) == len(perm)
        return tuple(pleft), tuple(pright)

    def compute_subgroups(self, only_full_groups=False):
        """
        Compute the subgroups of `self`.
        The option `only_full_groups` means that the left and right groups
        must be the left and right groups of the whole group.
        """
        self.ensure_pgroup()
        self.sgr0 = self.pgroup.subgroups()
        # print(len(self.sgr0),"subgroups in total")
        if self.achiral:
            self.sgr1 = self.sgr0
        else:
            self.sgr1 = [sgr for sgr in self.sgr0 if self.Pminus_1_minus_1 in sgr]
            # print(len(self.sgr1),"subgroups of them contain (-1,-1)")
        if only_full_groups:
            self.subgroups = [g for g in self.sgr1 if self.is_full_product(g)]
        else:
            self.subgroups = self.sgr1

    def is_full_product(self, g):
        if self.achiral:
            g.AxB = self.convert_Pgroup_to_AxB(g)
            L0 = {e.l for e in g.AxB}; L = L0.copy()
            L.update(self.Q1.minus_1 * e for e in L0)
            if len(L) != self.nL: return False
            R0 = {e.r for e in g.AxB}; R = R0.copy()
            R.update(self.Q1.minus_1 * e for e in R0)
            if len(R)!= self.nR: return False
            return True
        L = {left for e in g for left,right in (self.split(e),)}
        if len(L)!= self.nL: return False
        R = {right for e in g for left,right in (self.split(e),)}   
        if len(R)!= self.nR: return False
        return True

    def print_subgroup_lengths(self):
        lens = sorted({g.order() for g in self.subgroups})
        print("Lengths =",lens, "Number =", len(self.subgroups))
        print("Multiplicities =",
        list((l,len(list(x for x in self.subgroups if x.order()==l))) for l in lens))

    def element_as_pair(self, e):
        """
        Represent an element of the pgroup as a pair of left and right
        rotations (only for chiral; the chiral groups have a dictionary for
        this)
        """
        left, right = self.split(e)
        l = [x-1 for x in Permutation(left)]
        l += range(len(l),self.nL)
        r = [x-self.nL-1 for x in Permutation(right)[self.nL:]]             
        r += range(len(r),self.nR)           
        #print (l,r)
        lname = self.Q1.name_from_perm[tuple(l)]
        rname = self.Q2.name_from_perm[tuple(r)]
        return FDR(lname,rname)
    
    def convert_Pgroup_to_AxB(self, sgr):
        """
        Convert a permutation group `sgr` (as list of permutations)
        to a list of elements of self.
        """
        if self.achiral:
            g0 = [self.to_AxB_group[e.tuple()] for e in sgr]
            return g0 + [self.minus_1_minus_1*e for e in g0]
        else:
            self.Q1.ensure_name_from_perm()
            self.Q2.ensure_name_from_perm()
            return [self.element_as_pair(e) for e in sgr]

class AxB_product_group(AxB_group):
    """
    Construct AxB group whose left and right groups are `Q1` and `Q2`.
    """
    
    def __init__(self, Q1, Q2=None, name=None, achiral=False, limit=None):
        if achiral:
            assert Q2 is None
            Q2 = Q1
        self.Q1 = Q1
        self.Q2 = Q2
        self.achiral = achiral
        gens = [FDR(l,r)
                for l in [Q1.one]+list(Q1.gens)
                for r in [Q2.one]+list(Q2.gens)] [1:] # omit trivial generator (1,1)
        if achiral:
            gens.append(FDR(Q1.one, Q1.one, True))

        AxB_group.__init__(self, gens, name, limit)
        
    def find_equal_pairs(self, sgr):
        """
        Return a list of `p`'s such that `FDR(p, p)` in `sgr`.
        """
        return [p for p in self.Q1.q if FDR(p,p) in sgr.AxB]

    def find_kernels(self, sgr=None):
        """
        Set the left and the right kernel
        L0 = {p | [p, 1] in sgr}
        R0 = {p | [1, p] in sgr}
        as attributes to `sgr` (permutation group).
        If `sgr` is `None`, set them for `self`.
        """
        if self.achiral:
            raise NotImplementedError
            return
        if sgr is None:
            gr = self.pgroup
            who = self
        else:
            gr = who = sgr
            who 
        if not hasattr(gr, "AxB"):
            gr.AxB = self.convert_Pgroup_to_AxB(gr)
        who.L0 = list(p for p in self.Q1.q if FDR(p,self.Q2.one) in gr.AxB)
        who.R0 = list(q for q in self.Q2.q if FDR(self.Q1.one,q) in gr.AxB)

## Quaternion Groups
group_2T = Q_group([i, w])
group_2I = Q_group([i_I, w])
group_2O = Q_group([i_O, w])
group_2D2n = lambda n: Q_group([e(n), j])
group_2Cn = lambda n: Q_group([e(n)])

def extend_group(sgrp, gens, one=Q(1), limit=30000):
    ## NOT USED ANYMORE
    """
    Return the group generated by elements in `sgrp + gens`.
    `sgrp` is supposed be be already a subgroup, i.e. closed under "*".
    """
    def candidates(old, new):
        for A in new:
            for B in old:
                yield A*B
                yield B*B
            for B in new:
                yield A*B

    old_G = list(sgrp)
    # maintain group as a list because group 
    # elements are generated in a unique order
    G = list(sgrp) + list(gens) 
    # G (list) and G_set (set) always
    # contain the same elements
    G_set = set(G) 
    assert len(G)==len(G_set)
    new_in_G = gens
    while new_in_G:
        to_come = []
        for new_g in candidates(old_G,new_in_G):
            if new_g not in G_set:
                G.append(new_g)
                G_set.add(new_g)
                to_come.append(new_g)
                if len(G) > limit:
                    raise ValueError("group size over the limit")
        old_G.extend(new_in_G)
        new_in_G = to_come
    return G

class DisjointSet(object):
    # https://stackoverflow.com/a/3067672

    def __init__(self,elts=[]):
        self.leader = {a:a for a in elts} # maps a member to the cluster's leader
        self.cluster = {a:set([a]) for a in elts} # maps a cluster leader to the cluster (which is a set)

    def is_equivalent(self, a, b):
        leadera = self.leader.get(a)
        leadera = self.leader.get(b)
        return leadera is not None and leaderb is not None and leadera == leaderb
    
    def add(self, a, b):
        leadera = self.leader.get(a)
        leaderb = self.leader.get(b)
        if leadera is not None:
            if leaderb is not None:
                if leadera == leaderb: return # nothing to do
                clustera = self.cluster[leadera]
                clusterb = self.cluster[leaderb]
                if len(clustera) < len(clusterb):
                    a, leadera, clustera, b, leaderb, clusterb = b, leaderb, clusterb, a, leadera, clustera
                clustera |= clusterb
                del self.cluster[leaderb]
                for k in clusterb:
                    self.leader[k] = leadera
            else:
                self.cluster[leadera].add(b)
                self.leader[b] = leadera
        else:
            if leaderb is not None:
                self.cluster[leaderb].add(a)
                self.leader[a] = leaderb
            else: # works also for a==b
                self.leader[a] = self.leader[b] = a
                self.cluster[a] = set([a, b])

def cosets(gr, sgr):
    """
    The right cosets of subgroup `sgr` within `gr`.
    Here `gr` and `sgr` are sets or lists
    """
    def update_rem(sgr):
        for g in sgr:
            rem[g] = False
    rem = {g: True for g in gr}
    update_rem(sgr)
    cosets = [sgr]
    for u in list(gr):
        if not rem[u]: continue
        coset = [e*u for e in sgr]
        cosets.append(coset)
        update_rem(coset)
    return cosets

def time0():
    global s_time
    s_time = time()

def time1(command, file=sys.stdout):
    print ("Elapsed time:",time()-s_time,"seconds", flush=True, file=file)
    return command

###################################
## Functions to generate catalog ##
###################################
def CS_name_to_latex(name):
    """
    Turn CS `name` of a group to latex string.
    """
    def factor_to_latex(f):
        if f=="T'":
            return r"\overline T"
        if f=="O'":
            return r"\overline O"
        if f=="I'":
            return r"\overline I"
        if f.startswith("D'"): 
            return r"\overline D_{"+f[2:]+"}"  
        if f[0]in "CD": 
            power = ""
            if "^" in f:
                f, power = f.split("^")
                power = "^{("+power+")}"
            return f[0]+"_{"+f[1:]+"}" + power
        return f 

    assert name[0]=="+" 
    if name.startswith("+-"): 
        out = r"\pm" 
        name = name[2:] 
    else: 
        out = r"+" 
        name = name[1:] 

    frac, rest = name.split("[") 
    if frac != "": 
        assert frac.startswith("1/") 
        out += r"\frac1{"+frac[2:]+"}" 

    prod, rest = rest.split("]") 
    f1, f2 = prod.split("x") 
    out += "["+factor_to_latex(f1)+r"\times "+factor_to_latex(f2)+"]"
    if rest:
        assert rest.startswith(".")
        out += r"\cdot "
        if rest[1:]=="-2":
            out += r"\bar2"
        else:
            out += rest[1:]
    return out 

def quaternion_to_latex(q):
    """
    Turn a quaternion to latex string.
    """
    tr_table = {}
    for (e,l) in ((w, r'\omega'),(i_O, 'i_O'),(i_I, 'i_I'),(i_I2, "i'_I"),
                  (i_I_dag, r"i_I^\dag"), (i_I_dag*i_O*i, r"i_I^\dag i_O i"),              
                  (i_I*i_O*i, r"i_I i_O i"), (w*i,r'\omega i'),
                  (i*w,r'i\omega'), (w_bar, r'\overline \omega')):
        tr_table[e]=l
        tr_table[-e]="-"+l
    return tr_table.get(q,q)

def insert_into_catalog(gr, names_only=False, latex_only=False, verbose=False):
    """
    Insert group into `catalog`.
    The keys in the catalog are the fingerprints of the groups. That is,
    the output of the function `classify_group`.
    The values in the catalog are the corresponding AxB groups by default.
    """
    if verbose: print(f"processing {gr.name}. order: {gr.order//2:3}")

    value = gr
    if names_only: value = gr.name
    if latex_only: value = gr.latex
    clas = gr.classify_group()
    if clas in catalog:
        if value in catalog[clas]:
            print(f"WARNING: identical groups found:")
        else:
            print(f"WARNING: duplications found:")
        if names_only or latex_only:
            print([value] + [name for name in catalog[clas]])
        else:
            print([value.name] + [grp.name for grp in catalog[clas]])
        if value in catalog[clas]:
            print(f"WARNING: identical groups")
    catalog.setdefault(clas,[]).append(value)
    return 

def check_some_duval_groups():
    """
    Generate some of DuVal groups and compare them to the groups in `catalog`.
    Here `catalog` is assumed to be generated using `names_only`.
    """
    DuVal22 = mk_group([FDR(i,Q(1)),FDR(j,Q(1)),FDR(Q(1),i),FDR(Q(1),j),
                        FDR(w,w),FDR(w,-w)])
    DuVal41 = mk_group([FDR(i,Q(1)),FDR(j,Q(1)),FDR(Q(1),i),FDR(Q(1),j),
                        FDR(w,w),FDR(w,-w),FDR(Q(1),Q(1),1)])
    DuVal41b = mk_group([FDR(i,Q(1)),FDR(j,Q(1)),FDR(Q(1),i),FDR(Q(1),j),
                        FDR(w*w,w*w,1)])
    w2 = (1-i+j-k)/2
    DuVal42 = mk_group([FDR(i,Q(1)),FDR(j,Q(1)),FDR(Q(1),i),FDR(Q(1),j),
                        FDR(w2,w2.conjugate(),1), 
                        FDR(w,w.conjugate(),1)] # same coset
                        # one single reflection is not enough.
                      )
    for it, gr in [(22, DuVal22), (41, DuVal41), ("41b", DuVal41b),
                   ("42", DuVal42)]:
        print (f"DuVal{it!s:5} {len(gr):3}", catalog[classify_group(gr)])
    ## expected output:
    # DuVal22    192 ['+-1/3[TxT]']
    # DuVal41    384 ['+-1/3[TxT].2']
    # DuVal41b   384 ['+-1/3[TxT].2']
    # DuVal42    384 ["+-1/3[TxT'].2"]

    DuVal21 = mk_group([FDR(w,w),FDR(i,i),FDR(Q(1),Q(-1))])
    DuVal39 = mk_group([FDR(w,w),FDR(i,i),FDR(Q(1),Q(-1)),
                       FDR(Q(1),Q(1),1)])
    DuVal40 = mk_group([FDR(w,w),FDR(i,i),FDR(Q(1),Q(-1)),
                       FDR(i_O,i_O,1)])
    DuVal39p = mk_group([FDR(w,w),FDR(i,i), FDR(Q(1),Q(1),1)])
    DuVal40p = mk_group([FDR(w,w),FDR(i,i), FDR(i_O,i_O,1)])
    DuVal39pp = mk_group([FDR(w,w),FDR(i,i), FDR(w,-w,1)])
    DuVal40pp = mk_group([FDR(w,w),FDR(i,i), FDR(i_O,-i_O,1)])

    for it, gr in [(21, DuVal21), (39, DuVal39), (40, DuVal40),
                   ("39'", DuVal39p), ("40'", DuVal40p),
                   ("39'-", DuVal39pp), ("40'-", DuVal40pp)]:
        print ("DuVal",f"{it!s:5}", catalog[classify_group(gr)])
    ## expected output
    # DuVal 21    ['+-1/12[TxT]']
    # DuVal 39    ['+-1/12[TxT].2']
    # DuVal 40    ["+-1/12[TxT'].2"]
    # DuVal 39'   ['+1/12[TxT].2_3']
    # DuVal 40'   ["+1/12[TxT'].2_1"]
    # DuVal 39'-  ['+1/12[TxT].2_1']
    # DuVal 40'-  ["+1/12[TxT'].2_3"]
    return

def check_Cs_suffix():
    """
    For each group in `catalog`, list if it contains *[1, 1] (suffix 2_3) or
    *[-1, 1] (suffix 2_1). Here `catalog` is assumed to be generated using
    `names_only`.
    """
    for clas, names in catalog.items():
        typ = ""
        if any(k.startswith("*0:") for k in clas.split()):
            typ += "contains 2_3   "
        if any(k.startswith("*1:") for k in clas.split()):
            typ += "contains 2_1   "
        if typ == "":
            clas = "chiral"
        print(f"{' '.join([name for name in names]):20} {typ}")
    return

def generate_polyhedral_and_axial(names_only=False, latex_only=False,
                                  verbose=True):
    """
    Generate the polyhedral and the axial groups and insert them in `catalog`.

    INPUT:
    - `names_only` -- boolean: If True, inserted values in `catalogs` are
       names of the groups.
    - `latex_only` -- boolean: If True, inserted values in `catalogs` are 
       latex names of the groups. This option overrides `names_only`.
    - `verbose` -- boolean: More information.

    REMARK:
    `catalog` should be defined before running this function. Otherwise,
    you get a `NameError`.

    RUNNING TIME:
    46 groups in 181.6938283443451s.
    """
    ## check if catalog is already defined
    try:
        _ = len(catalog)
        if verbose: print(f"# catalog contains {len(catalog)} groups")
    except NameError:
        raise NameError("catalog is not defined")

    ##  measuring variables
    start_time = time()
    num_groups = 0

    ## generators for polyhedral 3-dimensional chiral groups
    enum_gr3s = ("I", [i_I, w]),("O", [i_O, w]),("T", [i, w])

    ## +-[LxR]
    for (gr3_1_name,gr3_1_gens) in enum_gr3s:
        for (gr3_2_name,gr3_2_gens) in enum_gr3s:
            L = mk_group(gr3_1_gens)
            R = mk_group(gr3_2_gens)
            group = [FDR(l,r) for l in L for r in R]
            name = f"+-[{gr3_1_name}x{gr3_2_name}]"
            gens = [FDR(l, Q(1)) for l in gr3_1_gens] + [FDR(Q(1), r) for r in gr3_2_gens]
            gr = AxB_group(gens=gens, name=name, group=group)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)
            num_groups += 1

    ## +-[LxL].2
    for (gr3_name, gr3_gens) in enum_gr3s:
        L = R = mk_group(gr3_gens)
        group = [FDR(l, r) for l in L for r in R]+[FDR(l, r, True) for l in L for r in R]
        name = f"+-[{gr3_name}x{gr3_name}].2"
        gens = [FDR(l, Q(1)) for l in gr3_gens] + [FDR(Q(1), r) for r in gr3_gens] + [FDR(Q(1), Q(1), True)]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)
        num_groups += 1

    ## +-1/60 IxI
    for ext_element, suffix_temp in [(FDR(Q(1), Q(1),0), ""),
                                     (FDR(Q(1), Q(1),1), ".2_3"),
                                     (FDR(Q(1), Q(-1),1), ".2_1")]:
       # print("reflection", suffix_temp, ", ".join((f"{unit}->{ext_element.act(unit)}" for unit in (1,i,j,k))))
        for sign, prefix in [(Q(-1), "+-"), (Q(1), "+")]:
            suffix = suffix_temp
            if prefix == "+-":
                if suffix_temp ==".2_1":
                    continue
                if suffix_temp ==".2_3":
                    suffix = ".2"
            for g, g_name in [(i_I, "I"), (i_I2, "I'")]:
                gens = [FDR(w,w), FDR(i_I, g), FDR(Q(1), sign), ext_element]
                name = f"{prefix}1/60[Ix{g_name}]{suffix}"
                gr = AxB_group(gens, name=name, limit=480)
                gr.latex = CS_name_to_latex(name)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1

    ## 1/24 OxO
    for ext_element, suffix_temp in [(FDR(Q(1), Q(1), 0), ""),
                                     (FDR(Q(1), Q(1), 1), ".2_3"),
                                     (FDR(Q(1), Q(-1), 1), ".2_1")]:
        for sign, prefix in [(Q(-1), "+-"), (Q(1), "+")]:
            suffix = suffix_temp
            if prefix == "+-":
                if suffix_temp == ".2_1":
                    continue
                elif suffix_temp == ".2_3":
                    suffix = ".2"
            for g, g_name in [(i_O, "O"), (-i_O, "O'")]:
                if g_name == "O'" and prefix == "+-":
                    continue
                gens = [FDR(w,w),FDR(i_O, g), FDR(Q(1), sign), ext_element]
                name = f"{prefix}1/24[Ox{g_name}]{suffix}"
                gr = AxB_group(gens, name=name, limit=480)
                gr.latex = CS_name_to_latex(name)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1

    ## 1/2 OxO and 1/6 OxO
    gensOxO = (
        ("1/2",
         [FDR(i, Q(1)), FDR(w, Q(1)), FDR(Q(1), i), FDR(Q(1), w), FDR(i_O, i_O)],
         [("", []), (".2", [FDR(Q(1), Q(1), 1)]), (".-2", [FDR(Q(1), i_O, 1)])]),
        ("1/6",
         [FDR(i, Q(1)), FDR(j, Q(1)), FDR(Q(1), i),
          FDR(Q(1), j), FDR(w, w), FDR(i_O, i_O)],
         [("", []), (".2", [FDR(Q(1), Q(1), 1)])])
    )
    for frac, gens0, exts in gensOxO:
        for suffix, ext_element in exts:
            name = "+-" + frac + "[OxO]"+suffix
            gr = AxB_group(gens0+ext_element, name=name, limit=3000)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)
            num_groups += 1

    ## 1/12 TxT
    for ext_element, suffix_temp in [(FDR(Q(1),Q(1),0), ""),
                                     (FDR(Q(1),Q(1),1), ".2_3"),
                                     (FDR(Q(1),Q(-1),1), ".2_1")]:
        for sign, prefix in [(Q(-1), "+-"), (Q(1), "+")]:
            suffix = suffix_temp
            if prefix =="+-":
                if suffix_temp == ".2_1":
                    continue
                elif suffix_temp ==".2_3":
                    suffix = ".2"
    #       for g1, g2, g_name in [(w,i,"T"),(w.conjugate(),-i,"T'")]:
            for g1, g2, g_name in [(w, i, "T"), (w, -i, "T'")]:
                if suffix_temp == "" and g_name == "T'":
                    continue                
                gens = [FDR(w, g1), FDR(i, g2), FDR(Q(1), sign), ext_element]            
                ## generators for TxT' different from CS:
                if g_name=="T'":
                    if prefix == "+-":
                        gens = [FDR(w, -w), FDR(i_O, i_O,1)]
                    elif suffix_temp ==".2_1":
                        gens = [FDR(w, w), FDR(i_O, i_O,1)]
                    elif suffix_temp ==".2_3":
                        gens = [FDR(w, w), FDR(i_O, -i_O,1)]
                    else:
                        raise ValueError
                name = f"{prefix}1/12[Tx{g_name}]{suffix}"
                gr = AxB_group(gens, name=name, limit=480)
                gr.latex = CS_name_to_latex(name)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1

    ## 1/3 TxT
    for ext_element, suffix in [(FDR(Q(1), Q(1), 0), ""),
                                (FDR(Q(1), Q(1), 1), ".2")
                               ]:
        for sign, prefix in [(Q(-1),"+-"),]:
            for g1, g2, g_name in [(w, i, "T"), (w.conjugate(), -i, "T'")]:
                if g_name == "T'" and suffix == "": continue
                gens = [FDR(i, Q(1)), FDR(j, Q(1)), FDR(Q(1), i), FDR(Q(1), j),
                        FDR(w, g1), FDR(Q(1), sign), ext_element]          
                name = f"{prefix}1/3[Tx{g_name}]{suffix}"
                gr = AxB_group(gens, name=name, limit=480)
                gr.latex = CS_name_to_latex(name)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1
    
    print("finished generating polyhedral and axial groups")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")
    return

def toroidal_group(category, p1=None, p2=None, p3=None):
    """
    Return the toroidal group with the given parameters according to the
    classification of these groups in the article "Towards a geometric
    understanding of 4-dimensional point groups."
    """
    assert category in "1./\\X|+L*", "unknown category"

    def T(phi1, phi2):
        """
        a torus translation T_{phi1, phi2} to a 4-dimensional rotation FDR
        """
        return FDR(quat_2D(-(1/2*phi1 + 1/2*phi2)), quat_2D(1/2*phi1 - 1/2*phi2))
    
    if category == "1":
        m, n, s = p1, p2, p3
        gens = [T(2/m, 2/m), T(2/n + 2*s/(m*n), 2*s/(m*n))]
        g = AxB_group(gens, name=f"G1^({s})_{m},{n}")
        g.latex_placeholder = r'\grp'+category+'_{{{1},{2}}}^{{({3})}}'
        if (n == 1) or (n == 2 and m % 2 == 1):
            g.latex_placeholder = r'\grp'+category+'_{{{1},{2}}}'
        assert g.order == 2*m*n

    elif category == ".":
        m, n, s = p1, p2, p3
        gens = [T(2/m, 2/m),
                T(2/n + 2*s/(m*n), 2*s/(m*n)),
                FDR(quat_j, quat_j)]
        g = AxB_group(gens, name=f"G.^({s})_{m},{n}")
        g.latex_placeholder = r'\grp'+category+'_{{{1},{2}}}^{{({3})}}'
        if (n == 1) or (n == 2 and m % 2 == 1):
            g.latex_placeholder = r'\grp'+category+'_{{{1},{2}}}'
        assert g.order == 4*m*n

    elif category in "\\/":
        wallpaper, m, n = p1, p2, p3
        gens = [FDR(quat_2D(2/m), quat_1), FDR(quat_1, quat_2D(2/n))]
        if category == "/":
            mirror = FDR(quat_i, quat_k)
        else:
            mirror = FDR(quat_k*quat_minus1, quat_i)
        if wallpaper == "pm":
            assert m % 2 == n % 2 == 0
            gens += [mirror]
            multiplier = 1
        elif wallpaper == "pg":
            assert m % 2 == n % 2 == 0
            if category == "/":
                gens += [FDR(quat_2D(1/m), quat_1)*mirror]
            else:
                gens += [FDR(quat_1, quat_2D(1/n))*mirror]
            multiplier = 1
        elif wallpaper == "cm":
            assert (m - n) % 2 == 0
            gens += [FDR(quat_2D(1/m), quat_2D(1/n)), mirror]
            multiplier = 2
        else:
            raise ValueError(wallpaper)
        g = AxB_group(gens, name=f"G{category}^{wallpaper}_{m},{n}")
        if category == "/":
            g.latex_placeholder = r'\grp/'
        else:
            g.latex_placeholder = r'\grp\setminus'
        g.latex_placeholder += r'^{{\mathbf{{'+wallpaper+r'}}}}_{{{2},{3}}}'
        assert g.order == 2*multiplier*m*n

    elif category == "X":
        wallpaper, m, n = p1, p2, p3
        gens = [FDR(quat_2D(2/m), quat_1), FDR(quat_1, quat_2D(2/n))]
        mirror1 = FDR(quat_i, quat_k)                # /
        mirror2 = FDR(quat_k*quat_minus1, quat_i)    # \
        if wallpaper == "cmm":
            assert (m - n) % 2 == 0
            multiplier = 4
        else:
            assert m % 2 == n % 2 == 0
            multiplier = 2
        if wallpaper == "cmm":
            gens += [FDR(quat_2D(1/m), quat_2D(1/n)), mirror1, mirror2]
        elif wallpaper == "pmm":
            gens += [mirror1, mirror2]
        elif wallpaper == "pgg":
            gens += [FDR(quat_2D(1/m), quat_2D(1/n))*mirror1,
                     FDR(quat_2D(1/m), quat_2D(1/n))*mirror2]
        elif wallpaper == "pmg":
            gens += [FDR(quat_1, quat_2D(1/n))*mirror1,
                     FDR(quat_1, quat_2D(1/n))*mirror2]
        elif wallpaper == "pgm":
            gens += [FDR(quat_2D(1/m), quat_1)*mirror1,
                     FDR(quat_2D(1/m), quat_1)*mirror2]
        else:
            raise ValueError(wallpaper)
        g = AxB_group(gens, name=f"G{category}^{wallpaper}_{m},{n}")
        wallpaper_full = wallpaper[:1]+"2"+wallpaper[1:]  # according to new IT
        g.latex_placeholder = r'\grp X^{{\mathbf{{'+wallpaper_full+r'}}}}_{{{2},{3}}}'
        assert g.order == 2*multiplier*m*n

    elif category == "|":
        wallpaper, m, n = p1, p2, p3
        gens = [T(-2/n, 0), T(0, -2/m)]
        h_reflection = FDR(quat_i, quat_i, True)            #|m
        if wallpaper == "pm":
            gens += [h_reflection]
            multiplier = 2
        elif wallpaper == "cm":
            gens += [T(-1/n, -1/m), h_reflection]
            multiplier = 4
        elif wallpaper == "pg":
            gens += [h_reflection*T(0, -1/m)]
            multiplier = 2
        else:
            raise ValueError(wallpaper)
        g = AxB_group(gens, name=f"G{category}^{wallpaper}_{m},{n}")
        g.latex_placeholder = r'\grp'+category+ r'^{{\mathbf{{'+wallpaper+r'}}}}_{{{2},{3}}}'
        assert g.order == 2*multiplier*m*n
    
    elif category == "+":
        wallpaper, m, n = p1, p2, p3
        gens = [T(-2/n, 0), T(0, -2/m)]
        h_reflection = FDR(quat_i, quat_i, True)            #|m
        v_reflection = FDR(quat_k, quat_k, True)            #-m
        if wallpaper == "pmm":
            gens += [h_reflection, v_reflection]
            multiplier = 4
        elif wallpaper == "pmg":    # mirror is vertical
            gens += [h_reflection*T(-1/n, 0), v_reflection*T(-1/n, 0)]
            multiplier = 4
        elif wallpaper == "pgg":
            gens += [h_reflection*T(-1/n, -1/m), v_reflection*T(-1/n, -1/m)]
            multiplier = 4
        elif wallpaper == "cmm":
            gens += [T(-1/n, -1/m), h_reflection, v_reflection]
            multiplier = 8
        g = AxB_group(gens, name=f"G{category}^{wallpaper}_{m},{n}")
        wallpaper_full = wallpaper[:1]+"2"+wallpaper[1:] # according to new IT
        g.latex_placeholder = r'\grp'+category+ r'^{{\mathbf{{'+wallpaper_full+r'}}}}_{{{2},{3}}}'
        assert g.order==2*multiplier*m*n

    elif category == "L":
        a, b = p1, p2
        assert p3 is None
        wallpaper = "p4"
        gens = [T(2*a/(a^2+b^2), 2*b/(a^2+b^2)),
                T(2*b/(a^2+b^2), -2*a/(a^2+b^2))]
        gens += [FDR(quat_minus1*quat_j, quat_1, True)]
        #g = AxB_group(gens, name=f"G{category}^{wallpaper}_{a},{b}")
        g = AxB_group(gens, name=f"G{category}_{a},{b}")
        assert g.order == 2*4*(a^2+b^2)
        g.latex_placeholder = r'\grp L_{{{1},{2}}}'

    elif category == "*":
        wallpaper, n = p1, p2
        assert p3 is None
        h_reflection = FDR(quat_i, quat_i, True)            #|
        swap = FDR(quat_i, quat_k)                          #/
        if wallpaper in ["p4mS", "p4gS"]:
            gens = [FDR(quat_2D(1/n), quat_1), FDR(quat_1, quat_2D(1/n))]
            if wallpaper == "p4mS":
                gens += [h_reflection, swap]
            else:
                gens += [h_reflection*T(0, -1/n), swap*T(0,-1/n)]
            multiplier = 16
        elif wallpaper in ["p4mU", "p4gU"]:
            gens = [T(-2/n, 0), T(0, -2/n)]
            if wallpaper == "p4mU":
                gens += [h_reflection, swap]
            else:
                # gens += [T(1/n, 1/n)*h_reflection, swap] #center 2f-rotation
                gens += [h_reflection*T(-1/n, -1/n), swap*T(-1/n, -1/n)] #center 4f-rotation
            multiplier = 8
        else:
            raise ValueError(wallpaper)
        g = AxB_group(gens, name=f"G{category}^{wallpaper[2:]}_{n}") # cut out "p4"
        assert g.order == 2*multiplier*n*n, (g.order,2*multiplier*n*n)
        wallpaper_full = (r'\mathbf{{'+wallpaper[:3]
                          +'m'     # according to new IT
                          +r"}}\textrm{{"+wallpaper[-1]+"}}")
        g.latex_placeholder = r'\grp'+category+ r'^{{'+wallpaper_full+r'}}_{{{2}}}'

    g.latex = g.latex_placeholder.format(0,p1,p2,p3)
    return g

def generate_chiral_toroidal(orderbound=130, generate_full=False,
                             names_only=False, latex_only=False,
                             verbose=True):
    """
    Generate chiral toroidal groups and insert them in `catalog`.

    INPUT:
    - `orderbound` -- int: Order bound on the generated group.
    - `generate_full` -- boolean: Whether to include duplications or not.
    - `names_only` -- boolean: If True, inserted values in `catalogs` are
       names of the groups.
    - `latex_only` -- boolean: If True, inserted values in `catalogs` are 
       latex names of the groups. This option overrides `names_only`.
    - `verbose` -- boolean: More information.

    REMARK:
    `catalog` should be defined before running this function. Otherwise,
    you get a `NameError`.

    RUNNING TIME:
    orderbound | generate_full | num_groups | time
    ---------- | ------------- | ---------- | ------------------
    130        | True          | 10190      | 45.391656160354614s
    130        | False         | 9666       | 43.985748291015625s
    500        | True          | 136149     | 2862.3144364356995s
    500        | False         | 134099     | 2677.123522758484s
    1000       | False         | 527307     | 22028.41566133499s
    """
    ## check if catalog is already defined
    try:
        _ = len(catalog)
        if verbose: print(f"# catalog contains {len(catalog)} groups")
    except NameError:
        raise NameError("catalog is not defined")

    ##  measuring variables
    start_time = time()
    num_groups = 0

    print("start generating chiral toroidal groups up to order", orderbound)
    
    ## grp. and grp1
    for typ in "1.":
        for m in range(1, 1+orderbound):
            for n in range(1, 1+orderbound):
                if m*n > orderbound or (typ=="." and 2*m*n>orderbound): break
                for s in range(-(m//2),1+(n-m)//2):
                    gr = toroidal_group(typ,m,n,s)
                    if typ=="." and m<=2 and n==1:
                        clas = classify_group(gr.group)
                        assert clas in catalog
                        if not generate_full:
                            continue 
                    insert_into_catalog(gr, names_only, latex_only, verbose)
                    num_groups += 1

    ## grp\ and grp/
    for m in range(1,1+orderbound//2):
        for n in range(1,1+orderbound//2):
            if 2*m*n>orderbound: break
            for typ in "\\/":
                for wallp in ("pm","pg","cm"):
                    if wallp != "cm":
                        if 4*m*n>orderbound: continue
                        gr = toroidal_group(typ,wallp,2*m,2*n)
                    if wallp == "cm":
                        if (m-n)%2 != 0: continue
                        gr = toroidal_group(typ,wallp,m,n)
                    if ((wallp in ("pm","cm") and (n==1 or m==1))
                        or (wallp=="pg" and ((typ=="/" and n==1) or (typ=="\\" and m==1)))
                        or (wallp=="cm" and ((typ=="/" and n==2) or (typ=="\\" and m==2)))
                        ): 
                        clas = classify_group(gr.group)
                        assert clas in catalog
                        if not generate_full:
                            continue 
                    insert_into_catalog(gr, names_only, latex_only, verbose)
                    num_groups += 1

    ## grpX
    for m in range(1,1+orderbound//4):
        for n in range(1,1+orderbound//4):
            if 4*m*n>orderbound: break
            for wallp in ("pmm","pgg","cmm","pmg","pgm"):
                    if wallp == "cmm":
                        if (m-n)%2 != 0: continue
                        gr = toroidal_group("X",wallp,m,n)
                    else:
                        if 8*m*n>orderbound: continue
                        gr = toroidal_group("X",wallp,2*m,2*n)
                    if n==1 or m==1 or (wallp=="cmm" and (n==2 or m==2)):                     
                        clas = classify_group(gr.group)
                        assert clas in catalog
                        if not generate_full:
                            continue 
                    insert_into_catalog(gr, names_only, latex_only, verbose)
                    num_groups += 1
    
    print("finished generating chiral toroidal groups")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")

    return

def generate_achiral_toroidal(orderbound=130, generate_full=False,
                              names_only=False, latex_only=False,
                              verbose=True):
    """
    Generate achiral toroidal groups and insert them in `catalog`.

    INPUT:
    - `orderbound` -- int: Order bound on the generated group.
    - `generate_full` -- boolean: Whether to include duplications or not.
    - `names_only` -- boolean: If True, inserted values in `catalogs` are
       names of the groups.
    - `latex_only` -- boolean: If True, inserted values in `catalogs` are 
       latex names of the groups. This option overrides `names_only`.
    - `verbose` -- boolean: More information.

    REMARK:
    `catalog` should be defined before running this function. Otherwise,
    you get a `NameError`.

    RUNNING TIME:
    orderbound | generate_full | num_groups | time
    ---------- | ------------- | ---------- | ------------------
    130        | True          | 1131       | 5.745510816574097s
    130        | False         | 973        | 4.516056537628174s
    500        | True          | 5725       | 126.6421959400177s
    500        | False         | 4929       | 107.031977891922s
    1000       | False         | 11110      | 443.3721535205841s
    """
    ## check if catalog is already defined
    try:
        _ = len(catalog)
        if verbose: print(f"# catalog contains {len(catalog)} groups")
    except NameError:
        raise NameError("catalog is not defined")

    ##  measuring variables
    start_time = time()
    num_groups = 0

    print("start generating achiral toroidal groups up to order", orderbound)

    ## grp|
    for wallp in ("pm", "cm", "pg"):
        for m in IntegerRange(1, orderbound+1):
            for n in IntegerRange(1,orderbound+1):
                if wallp == "cm" and 4*m*n > orderbound: break
                else:
                    if 2*m*n > orderbound: break
                gr = toroidal_group("|", wallp, m, n)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1
   
    ## grp+
    for wallp in ("pmm", "pmg", "pgg", "cmm"):
        for m in IntegerRange(1, orderbound+1):
            for n in IntegerRange(1, orderbound+1):
                if wallp == "cmm" and 8*m*n > orderbound: break
                else:
                    if 4*m*n > orderbound: break
                if not generate_full:
                    if (m == 1 and n == 1) or ((m < n) and (wallp !="pmg")):
                        continue
                gr = toroidal_group("+", wallp, m, n)
                insert_into_catalog(gr, names_only, latex_only, verbose)
                num_groups += 1
   
    ## grpL
    for a in IntegerRange(1,orderbound+1):
        for b in IntegerRange(0,orderbound+1):
            if 4*(a^2+b^2) > orderbound: break
            if not generate_full:
                if ((a,b) in [(1,0), (1,1), (2,0)]) or (a < b):
                    continue
            gr = toroidal_group("L", a, b)
            insert_into_catalog(gr, names_only, latex_only, verbose)
            num_groups += 1
   
    ## grp*
    for wallp in ("p4mU", "p4gU", "p4mS", "p4gS"):
        for n in IntegerRange(1,orderbound+1):
            if wallp in ("p4mS", "p4gS") and 16*n*n > orderbound: break
            else:
                if 8*n*n > orderbound: break
            if not generate_full:
                if (wallp in ("p4mS", "p4gS") and n == 1) or (wallp in ("p4mU", "p4gU") and (n == 1 or n == 2)):
                    continue
            gr = toroidal_group("*", wallp, n)
            insert_into_catalog(gr, names_only, latex_only, verbose)
            num_groups += 1
    
    print("finished generating achiral toroidal groups")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")

    return

def generate_tubical(orderbound=1000, generate_full=False,
                     names_only=False, latex_only=False,
                     verbose=True):
    """
    Generate tubical groups and insert them in `catalog`.

    INPUT:
    - `orderbound` -- int: Order bound on the generated group.
    - `generate_full` -- boolean: Whether to include duplications or not.
    - `names_only` -- boolean: If True, inserted values in `catalogs` are
       names of the groups.
    - `latex_only` -- boolean: If True, inserted values in `catalogs` are 
       latex names of the groups. This option overrides `names_only`.
    - `verbose` -- boolean: More information.

    REMARK:
    `catalog` should be defined before running this function. Otherwise,
    you get a `NameError`.

    RUNNING TIME:
    orderbound | num_groups | time
    1000       | 418        | 337.77567744255066s
    2000       | 860        | 1966.2435729503632s
    5000       |            |
    """
    ## check if catalog is already defined
    try:
        _ = len(catalog)
        if verbose: print(f"# catalog contains {len(catalog)} groups")
    except NameError:
        raise NameError("catalog is not defined")

    ##  measuring variables
    start_time = time()
    num_groups = 0

    ## There are only 5 duplications among tubical groups:
    ## +-[C2xI] .= +-[C1xI] and +-[IxC2] .= +-[IxC1]
    ## +-[C2xO] .= +-[C1xO] and +-[OxC2] .= +-[OxC1]
    ## +-[C2xT] .= +-[C1xT] and +-[TxC2] .= +-[TxC1]
    ## +-1/2[OxD2] .= +-1/2[OxC2] and +-1/2[D2xO] .= +-1/2[C2xO]
    ## +-1/2[OxD4] .= +-1/2[OxD'4] and +-1/2[D4xO] .= +-1/2[D'4xO]
    ## to element them, it is enough to start n at some cases from 2
    starting_n = 1 if generate_full else 2

    print("start generating tubical groups up to order", orderbound)

    for (gr3_name, gr3_gens) in [("I", [i_I, w]), ("O", [i_O, w]), ("T",[i, w])]:
        gr3 = mk_group(gr3_gens)
        gr3_order = len(gr3)//2 
        ## +-[LxCn] and +-[CnxR]
        for n in range(1, orderbound//(2*gr3_order)+1): # 2 because of +-
            cyclic = Cn_group(n)
            
            group = [FDR(l, r) for l in cyclic.q for r in gr3]
            gens = [FDR(quat_1, r) for r in gr3_gens] + [FDR(quat_2D(1/n), Q(1))]
            name = f"+-[C{n}x{gr3_name}]"
            gr = AxB_group(gens=gens, name=name, group=group)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)

            group = [FDR(l, r) for r in cyclic.q for l in gr3]
            gens = [FDR(l, quat_1) for l in gr3_gens] + [FDR(Q(1), quat_2D(1/n))]
            name = f"+-[{gr3_name}xC{n}]"
            gr = AxB_group(gens=gens, name=name, group=group)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)

            num_groups += 2

        ## +-[LxD2n] and +-[D2nxR]
        for n in range(starting_n, orderbound//(4*gr3_order)+1):
            dihedral = Dn_group(2*n)
       
            group = [FDR(l, r) for l in dihedral.q for r in gr3]
            gens = [FDR(quat_1, r) for r in gr3_gens]
            gens += [FDR(quat_2D(1/n), Q(1)), FDR(quat_j, Q(1))]
            name = f"+-[D{2*n}x{gr3_name}]"
            gr = AxB_group(gens=gens, name=name, group=group)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)
        
            name = f"+-[{gr3_name}xD{2*n}]"
            gens = [FDR(e.r, e.l) for e in gens]
            group = [FDR(e.r, e.l) for e in group]
            gr = AxB_group(gens=gens, name=name, group=group)
            gr.latex = CS_name_to_latex(name)
            insert_into_catalog(gr, names_only, latex_only, verbose)
            
            num_groups += 2

    ## +-1/2[OxD2n] and +-1/2[D2nxO]
    factor = 48
    for n in range(starting_n, orderbound//factor+1): #n=1 : D2 becomes C2
        gens = [FDR(i, quat_1), FDR(w, quat_1), FDR(Q(1), quat_2D(1/n)),
                FDR(i_O, quat_j)]
        name1, name2 = "O", f"D{2*n}"
        name = f"+-1/2[{name1}x{name2}]"
        group = mk_group(gens, one=one_HQ, limit=4*24*n)
        assert len(group) == 2*factor*n
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        name = f"+-1/2[{name2}x{name1}]"
        gens = [FDR(e.r, e.l) for e in gens]
        group = [FDR(e.r, e.l) for e in group]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        num_groups += 2

    ## +-1/2[OxD'4n] and +-/12[D'4nxO]
    factor = 96        
    for n in range(starting_n, orderbound//factor+1): #n=1: Dbar4 becomes D4
        gens = [FDR(i, quat_1), FDR(w, quat_1), FDR(Q(1), quat_2D(1/n)),
                FDR(Q(1),quat_j), FDR(i_O, quat_2D(1/(2*n)))]
        name1, name2 = "O", f"D'{4*n}"
        name = f"+-1/2[{name1}x{name2}]"
        group = mk_group(gens, one=one_HQ, limit=8*24*n)
        assert(len(group) == 2*factor*n)
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        name = f"+-1/2[{name2}x{name1}]"
        gens = [FDR(e.r, e.l) for e in gens]
        group = [FDR(e.r, e.l) for e in group]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)
        
        num_groups += 2

    ## +-1/6[OxD6n] and +-1/6[D6nxO]
    factor = 48
    for n in range(1, orderbound//factor+1):
        gens = [FDR(i,quat_1),FDR(j,quat_1),FDR(Q(1),quat_2D(1/n)),
                FDR(i_O,quat_j), FDR(w, quat_2D(1/(3*n)))]
        name1, name2 = "O", f"D{6*n}"
        name = f"+-1/6[{name1}x{name2}]"
        group = mk_group(gens, one=one_HQ, limit=8*24*n)
        assert(len(group) == 2*factor*n)
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        name = f"+-1/6[{name2}x{name1}]"
        gens = [FDR(e.r, e.l) for e in gens]
        group = [FDR(e.r, e.l) for e in group]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)
        
        num_groups += 2
        
 
    ## +-1/2[OxC2n] and +-1/2[C2nxO]
    factor = 48
    for n in range(1, orderbound//factor+1):
        gens = [FDR(i, quat_1), FDR(w, quat_1), FDR(Q(1), quat_2D(1/n)),
                FDR(i_O, quat_2D(1/(2*n)))]
        name1, name2 = "O", f"C{2*n}"
        name = f"+-1/2[{name1}x{name2}]"
        group = mk_group(gens, one=one_HQ, limit=8*24*n)
        assert(len(group) == 2*factor*n)
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        name = f"+-1/2[{name2}x{name1}]"
        gens = [FDR(e.r, e.l) for e in gens]
        group = [FDR(e.r, e.l) for e in group]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)
        
        num_groups += 2
          
    factor = 24
    for n in range(1, orderbound//factor+1):
        gens = [FDR(i, quat_1), FDR(Q(1), quat_2D(1/n)),
                FDR(w, quat_2D(1/(3*n)))]
        name1, name2 = "T", f"C{3*n}"
        name = f"+-1/3[{name1}x{name2}]"
        group = mk_group(gens, one=one_HQ, limit=8*24*n)
        assert(len(group) == 2*factor*n)
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)

        name = f"+-1/3[{name2}x{name1}]"
        gens = [FDR(e.r, e.l) for e in gens]
        group = [FDR(e.r, e.l) for e in group]
        gr = AxB_group(gens=gens, name=name, group=group)
        gr.latex = CS_name_to_latex(name)
        insert_into_catalog(gr, names_only, latex_only, verbose)
        
        num_groups += 2

    print("finished generating tubical groups")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")

    return

print("done with loading of main_functions.sage")
