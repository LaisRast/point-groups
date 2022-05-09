#!/usr/bin/env python 
# coding: utf-8

## Checks if we have the correct conjugations for toroidal groups and
## generates a LaTeX table `tables/toroidal_duplications.tex` of them.

load("main_functions.sage")
load("conjugacy_Dn.sage")
from itertools import zip_longest

## timing
start_overall = time()

## conjugacy_Dn.sage overwrites e
e = lambda n: Q([cos(pi/n), sin(pi/n), 0, 0])

## testing bound on the parameters of groups with parameters
x_bound = 10
# x_bound = 20 # 2 minutes
# x_bound = 100 # 71115.83813118935s

## data to be checked and converted as a table
# Last entry of a tuple is None iff no parameter. Otherwise conditions.
data = [
    "LEFT",
        r"\hline",
        r"\multicolumn6{|l|}{chiral groups}\\",
        r"\hline",
        r"$\grp{1}^{{(s)}}_{m,n}$& $\grp{1}^{(s+n)}_{m,n}$& $[1,1]$ (equal) &"+
        r"$\grp{.}^{{(s)}}_{m,n}$& $\grp{.}^{(s+n)}_{m,n}$& $[1,1]$ (equal) \\",
        r"$\grp{1}^{{(s)}}_{m,n}$& $\grp{1}^{(-m-s)}_{m,n}$&$[i,k] = \sym{/}$ &"+
        r"$\grp{.}^{{(s)}}_{m,n}$& $\grp{.}^{(-m-s)}_{m,n}$&$[i,k] = \sym{/}$ \\",
    "RIGHT",
        (
            lambda m: (".", 1, 1, 0),
            lambda m: ("1", 1, 2, 0),
            lambda m: FDR(i + j, 1 + k),
            lambda m: None
        ),
        (
            lambda x: (".", 2, 1, -1),
            lambda x: ("1", 2, 2, 0),
            lambda x: FDR(i + j, i + j),
            lambda x: None
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("\\", "pm", 4*m-2, 2), 
            lambda m: (".", 4*m-2, 1, -2*m+1),
            lambda m: FDR(j + k, i + j),
            lambda m: m >= 1
        ),
        (
            lambda m: ("\\", "pm", 4*m, 2), 
            lambda m: (".", 4*m, 1, -2*m), 
            lambda m: FDR(1, i + j),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("\\", "pm", 2, 4*m-2),
            lambda m: ("1", 2, 4*m-2, 2*m-2),
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("\\", "pm", 2, 4*m),
            lambda m: ("1", 4, 2*m, -2),
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1
        ),
    "RIGHT",
        (
            lambda m: ("/", "pm", 2, 4*m-2),
            lambda m: (".", 2, 2*m-1, -1), 
            lambda m: FDR(i + j, j + k),
            lambda m: m >= 1
        ),
        (
            lambda m: ("/", "pm", 2, 4*m),
            lambda m: (".", 2, 2*m, -1), 
            lambda m: FDR(i + j, 1),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("/", "pm", 2*m, 2), 
            lambda m: ("1", 2*m, 2, 0),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("\\", "pg", 2, 4*m-2),
            lambda m: ("1", 4, 2*m-1, -2), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("\\", "pg", 2, 4*m), 
            lambda m: ("1", 2, 4*m, 2*m-1), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
    "RIGHT",
        (
            lambda m: ("/", "pg", 2*m, 2),
            lambda m: ("1", 2*m, 2, 1),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1
        ),
    "OUT",
        r"\hline",
        (   
            lambda m: ("\\", "cm", 2*m+1, 1),
            lambda m: (".", 2*m+1, 1, 0), 
            lambda m: FDR(j + k, 1 - k),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("\\", "cm", 1, 4*m-3), 
            lambda m: ("1", 1, 8*m-6, 2*m-2), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("\\", "cm", 1, 4*m-1), 
            lambda m: ("1", 1, 8*m-2, 2*m-1), 
            lambda m: FDR(1 - j, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("\\", "cm", 2, 4*m-2), 
            lambda m: ("\\", "pm", 4, 4*m-2), 
            lambda m: FDR(i + j, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("\\", "cm", 2, 4*m),
            lambda m: ("\\", "pg", 4, 4*m),
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
    "RIGHT",
        (
            lambda m: ("/", "cm", 1, 2*m+1),
            lambda m: (".", 1, 2*m+1, m),
            lambda m: FDR(i + j, j + k),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("/", "cm", 4*m-3, 1),
            lambda m: ("1", 4*m-3, 2, -2*m-2),
            lambda m: FDR(1, 1 - j),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("/", "cm", 4*m-1, 1),
            lambda m: ("1", 4*m-1, 2, -2*m+1),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("/", "cm", 4*m-2, 2), 
            lambda m: ("/", "pm", 4*m-2, 4), 
            lambda m: FDR(1, i + j),
            lambda m: m >= 1
        ),
        (
            lambda m: ("/", "cm", 4*m, 2),
            lambda m: ("/", "pg", 4*m, 4),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1 
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("X", "pmm", 2*m, 2),
            lambda m: (".", 2*m, 2, 0), 
            lambda m: FDR(1, i + k),
            lambda m: m >= 1 
        ),
    "RIGHT",
        (
            lambda m: ("X", "pmm", 2, 4*m-2), 
            lambda m: (".", 2, 4*m-2, 2*m-2),
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "pmm", 2, 4*m),
            lambda m: (".", 4, 2*m, -2), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("X", "pgm", 2*m, 2),
            lambda m: (".", 2*m, 2, 1),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1
        ),
	"EMPTY",
#        (
#            lambda m: ("X", "pgm", 4*m, 2),
#            lambda m: (".", 4*m, 2, 1),
#            lambda m: FDR(1, i + k), 
#            lambda m: m >= 1        # merge with above
#        ),
        (
            lambda m: ("X", "pgm", 2, 4*m-2),
            lambda m: ("/", "cm", 2, 4*m-2), 
            lambda m: FDR(i + j, j + k),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "pgm", 2, 4*m),
            lambda m: ("/", "pm", 4, 4*m), 
            lambda m: FDR(i + j, 1),
            lambda m: m >= 1 
        ),
    "RIGHT",
        (
            lambda m: ("X", "pmg", 2, 4*m-2),
            lambda m: (".", 4, 2*m-1, -2), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "pmg", 2, 4*m), 
            lambda m: (".", 2, 4*m, 2*m-1), 
            lambda m: FDR(i + k, 1), 
            lambda m: m >= 1 
        ),
        (
            lambda m: ("X", "pmg", 4*m-2, 2),
            lambda m: ("\\", "cm", 4*m-2, 2),
            lambda m: FDR(j + k, i + j),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "pmg", 4*m, 2), 
            lambda m: ("\\", "pm", 4*m, 4), 
            lambda m: FDR(1, i + j),
            lambda m: m >= 1 
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("X", "pgg", 4*m-2, 2), 
            lambda m: ("\\", "pm", 4*m-2, 4), 
            lambda m: FDR(j + k, i + j),
            lambda m: m >= 1
        ),
    "RIGHT",
        (
            lambda m: ("X", "pgg", 4*m, 2), 
            lambda m: ("\\", "cm", 4*m, 2),
            lambda m: FDR(1, i + j),
            lambda m: m >= 1 
        ),
    "LEFT",
        (
            lambda m: ("X", "pgg", 2, 4*m-2), 
            lambda m: ("/", "pm", 4, 4*m-2), 
            lambda m: FDR(i + j, j + k),
            lambda m: m >= 1
        ),
    "RIGHT",
        (
            lambda m: ("X", "pgg", 2, 4*m), 
            lambda m: ("/", "cm", 2, 4*m), 
            lambda m: FDR(i + j, 1),
            lambda m: m >= 1 
        ),
    "OUT",
        r"\hline",
        (
            lambda m: ("X", "cmm", 4*m-3, 1), 
            lambda m: (".", 4*m-3, 2, -2*m + 2),
            lambda m: FDR(1, 1 - j),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("X", "cmm", 4*m-1, 1),
            lambda m: (".", 4*m-1, 2, -2*m+1),
            lambda m: FDR(1, 1 + j),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("X", "cmm", 1, 4*m-3), 
            lambda m: (".", 1, 8*m-6, 2*m-2),
            lambda m: FDR(1 + j, 1),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "cmm", 1, 4*m-1), 
            lambda m: (".", 1, 8*m-2, 2*m-1),
            lambda m: FDR(1 - j, 1),
            lambda m: m >= 1 
        ),
    "RIGHT",
        (
            lambda m: ("X", "cmm", 4*m-2, 2), 
            lambda m: ("X", "pmm", 4*m-2, 4), 
            lambda m: FDR(j + k, i + j),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "cmm", 4*m, 2),
            lambda m: ("X", "pgm", 4*m, 4),
            lambda m: FDR(1, i + k),
            lambda m: m >= 1 
        ),
        (
            lambda m: ("X", "cmm", 2, 4*m-2), 
            lambda m: ("X", "pmm", 4, 4*m-2), 
            lambda m: FDR(i + j, j + k),
            lambda m: m >= 1
        ),
        (
            lambda m: ("X", "cmm", 2, 4*m), 
            lambda m: ("X", "pmg", 4, 4*m), 
            lambda m: FDR(i + k, 1),
            lambda m: m >= 1 
        ),
    "OUT",
        r"\hline",
        r"\multicolumn6{|l|}{achiral groups}\\",
        r"\hline",
        (
            lambda m,n: ("+", "pmm", m, n),
            lambda m,n: ("+", "pmm", n, m),
            lambda m,n: FDR(i, k), # or FDR(j, 1),
            lambda m,n: n > m
        ),
        (
            lambda m,n: ("+", "pgg", m, n),
            lambda m,n: ("+", "pgg", n, m),
            lambda m,n: FDR(i, k), # or FDR(j, 1),
            lambda m,n: n > m
        ),
        (
            lambda m,n: ("+", "cmm", m, n),
            lambda m,n: ("+", "cmm", n, m),
            lambda m,n: FDR(i, k), # or FDR(j, 1),
            lambda m,n: n > m
        ),
        (
            lambda m: ("+", "pmm", 1, 1),
            lambda m: ("|", "pm", 1, 2),
            lambda m: FDR(1 + k, 1 - k),
            lambda m: None
        ),
        (
            lambda m: ("+", "pmg", 1, 1),
            lambda m: ("|", "pm", 2, 1),
            lambda m: FDR(1 + k, i - j),
            lambda m: None
        ),
         (
            lambda m: ("+", "pgg", 1, 1),
            lambda m: ("|", "pg", 1, 2),
            lambda m: FDR(1 + k, 1 - k),
            lambda m : None
        ),
        (
            lambda m: ("+", "cmm", 1, 1),
            lambda m: ("|", "pm", 2, 2),
            lambda m: FDR(1 + k, 1 - k),
            lambda m: None
        ),
    "RIGHT",
        (
            lambda m,n: ("L", m, n),
            lambda m,n: ("L", n, m),
            lambda m,n: FDR(i,k), 
            lambda m,n: n > m
        ),
        (
            lambda m: ("L", 1, 0),
            lambda m: ("|", "pg", 2, 1),
            lambda m: FDR(1 + k, 1 - i + j + k),
            lambda m: None
        ),
        (
            lambda m: ("L", 1, 1),
            lambda m: ("|", "pg", 2, 2),
            lambda m: FDR(1 + k, 1 + i - j + k),
            lambda m: None
        ),
        (
            lambda m: ("L", 2, 0),
            lambda m: ("+", "pgg", 2, 2),
            lambda m: FDR(1 + k, 1 + k),
            lambda m: None
        ),
    "OUT",
        r"\hline",
    "LEFT",
        (
            lambda m: ("*", "p4mU", 1),
            lambda m: ("+", "pmg", 1, 2),
            lambda m: FDR(1 + k, 1 - i - j - k),
            lambda m: None
        ),
    "LEFT",
        (
            lambda m: ("*", "p4mU", 2),
            lambda m: ("+", "cmm", 2, 2),
            lambda m: FDR(1 + k, 1 + k),
            lambda m: None
        ),
    "RIGHT",
        (
            lambda m: ("*", "p4gU", 1),
            lambda m: ("+", "pgg", 2, 1),
            lambda m: FDR(1 + k, 1 + i - j + k),
            lambda m: None
        ),
    "RIGHT",
        (
            lambda m: ("*", "p4gU", 2),
            lambda m: ("L", 2, 2),
            lambda m: FDR(1 + j, 1 + j),
            lambda m: None
        ),
    "LEFT",
        (
            lambda m: ("*", "p4mS", 1),
            lambda m: ("+", "pmg", 2, 2),
            lambda m: FDR(1 + k, 1 + i + j - k),
            lambda m: None
        ),
    "RIGHT",
        (
            lambda m: ("*", "p4gS", 1),
            lambda m: ("|", "cm", 2, 2),
            lambda m: FDR(1 + k, 1 + k),
            lambda m: None
        ),
    "OUT",
]

## dictionary to convert common quaternions to quat_2D
## the keys are not normalized to match the conjugations
## given in data
to_quat_2D = {e.to_exact_quaternion(): e
              for r in IntegerRange(4) for p in (0,1)
              for e in (quat_2D(2*r/4, p),)}
to_quat_2D[1+i] = quat_2D(1/4)
to_quat_2D[j+k] = quat_2D(1/4, 1)
to_quat_2D[-1+i] = quat_2D(3/4)
to_quat_2D[-j+k] = quat_2D(3/4, 1)
to_quat_2D[-1-i] = quat_2D(5/4)
to_quat_2D[-j-k] = quat_2D(5/4, 1)
to_quat_2D[1-i] = quat_2D(7/4)
to_quat_2D[j-k] = quat_2D(7/4, 1)

## test to_quat_2D
for quaternion, its_quat_2D in to_quat_2D.items():
    assert quaternion/vector(quaternion).norm() == its_quat_2D.to_exact_quaternion()

def make_conjname(conj, OK):
    if (type(conj) == str and conj == "MIRROR") or \
       (type(conj) is FDR and conj == FDR(i,k)):
        conjname = str(FDR(i,k))+r"=\sym/"
    else:
        assert type(conj) == FDR
        conjname =  str(conj)
    if not OK:
        conjname += "WRONG"
    return conjname

def check_conj(gr1, gr2, conj, achiral, name):
    if type(conj) == str and conj == "MIRROR":
        conj = FDR(i, k)

    gr1_original, gr2_original = gr1.copy(), gr2.copy()
    conj_original= conj

    ## try to convert conj componenets to quat_2D
    ## using to_quat_2D dictionary.
    conj_l = to_quat_2D[conj.l] if conj.l in to_quat_2D else conj.l
    conj_r = to_quat_2D[conj.r] if conj.r in to_quat_2D else conj.r
    conj = FDR(conj_l, conj_r)

    try:
        if achiral:
            if type(conj.l) is quat_2D and type(conj.r) is quat_2D:
                gr1 = [AxB_group.conjugate_element(e, conj.l, conj.r) for e in gr1]
            else:
                for tr, con_alt in transform_D4xD4:
                    if con_alt[:2] in ((conj.l, conj.r),(-conj.l, -conj.r)):
                        gr1 = [tr[e] for e in gr1]
                        break
                else:
                    raise ValueError (conj, "not known")
        else:
            if type(conj.l) is not quat_2D:
                if conj.l != 1:
                   for tr, con_alt in transform_D4:
                       if con_alt[0] in (conj.l, -conj.l):
                           gr1 = [FDR(tr[e.l], e.r) for e in gr1]
                           break
                   else:
                        raise ValueError (conj.l, "not known")
            else:
                gr1 = [AxB_group.conjugate_element(e, conj.l) for e in gr1]

            if type(conj.r) is not quat_2D:
                if conj.r != 1:
                   for tr, con_alt in transform_D4:
                       if con_alt[0] in (conj.r, -conj.r):
                           gr1 = [FDR(e.l, tr[e.r]) for e in gr1]
                           break
                   else:
                        raise ValueError (conj.r, "not known")
            else:
                gr1 = [AxB_group.conjugate_element(e, quat_1, conj.r) for e in gr1]
    except (ValueError,KeyError):
        print("fall back to exact quaternions in", name, "for conj", conj)
        gr1 = [e.to_exact_quaternions() for e in gr1_original]
        gr2 = [e.to_exact_quaternions() for e in gr2_original]
        gr1 = [AxB_group.conjugate_element(e, conj_original.l, conj_original.r,
                                           1/Q(conj_original.l).reduced_norm(),
                                           1/Q(conj_original.r).reduced_norm()) for e in gr1]
    
    gr1_set = set(gr1)
    gr2_set = set(gr2)
    if gr1_set != gr2_set:
        print("WRONG", name, conj)
        # raise
        return False
    assert gr1_set == gr2_set, name
    return True

def process_entry(entry, it=0):
    """
    process the entries of data when entry is a 4-tuple
    """
    Gr1, Gr2, Conj, cond = entry

    ## figure out number of parameters
    try:
        if cond(1) is None:
            numparams = 0
        else:
            numparams = 1
    except TypeError:
        numparams = 2

    ## figure out parameters range
    if numparams == 0:
        par_range = [(1,)]
    elif numparams == 1:
        par_range = [(x,) for x in IntegerRange(1,x_bound)]
    elif numparams == 2:
        par_range = [(x,y) for x in IntegerRange(1,x_bound)
                           for y in IntegerRange(1,x_bound)]

    ## check conjugation
    for pars in par_range:
        if numparams > 0 and not cond(*pars): continue
        params1 = Gr1(*pars)
        params2 = Gr2(*pars)
        conj = Conj(*pars)

        gr1 = toroidal_group(*params1)
        gr2 = toroidal_group(*params2)
        gr1name = gr1.latex_placeholder.format(*(latex(x) for x in params1))
        gr2name = gr2.latex_placeholder.format(*(latex(x) for x in params2))
        assert gr1name == gr1.latex
        assert gr2name == gr2.latex

        ## attention: order of groups is gr2, gr1
        OK = check_conj(gr2.group, gr1.group, conj, gr1.achiral,
                        (it, gr1.order, gr1.name, gr2.name))
        conjname = make_conjname(conj, OK)

    ## create gr1 and gr2 names
    if numparams == 1:
        var("m")
        params1 = Gr1(m)
        params2 = Gr2(m)
        gr1name = gr1.latex_placeholder.format(*(latex(x) for x in params1))
        gr2name = gr2.latex_placeholder.format(*(latex(x) for x in params2))
        conjname = make_conjname(Conj(m), OK)
    if numparams == 2:
        var("m n a b")
        params1 = Gr1(m,n)
        params2 = Gr2(m,n)
        if params1[0]=="L":
            params1 = Gr1(a,b)
            params2 = Gr2(a,b)
        gr1name = gr1.latex_placeholder.format(*(latex(x) for x in params1))
        gr2name = gr2.latex_placeholder.format(*(latex(x) for x in params2))
        conjname = make_conjname(Conj(m,n),OK)
    return fr'${gr1name:25}$& ${gr2name:25}$& ${conjname}$'

## output dir
output_dir = "tables/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

## start writing
outfile = open(output_dir + "tbl_toroidal_duplications.tex", "w")
outfile.write(r"""%{
\begin{tabular}[htb]{|l|l|l||l|l|l|}
\hline
$G_1$ & $G_2$ & $[\hat{l}, \hat{r}]$%: $G_1=h^{-1}G_2h$
& $G_1$ & $G_2$& $[\hat{l}, \hat{r}]$%: $G_1=h^{-1}G_2h$
\\
""")

fix1 = '$[1 + j,1 + i + j - k]$'
fix2 = '$[1 + j,1 {+} i {+} j {-} k]$'

column = left_column = []
right_column = []
for it, entry in enumerate(data):
    if entry == "LEFT":
        column = left_column
    elif entry == "RIGHT":
        column = right_column
    elif entry == "EMPTY":
        column.append("&&")
    elif entry == "OUT":
        for l, r in zip_longest(left_column, right_column, fillvalue="&&"):
            print(l+"&"+r+r'\\', file=outfile)
        column = left_column = []
        right_column = []
    elif type(entry) is str:
        print(entry, file=outfile)
    else:
        line = process_entry(entry, it)
        if line.endswith(fix1):
            line = line[:-len(fix1)] + fix2
        column.append(line)
    
outfile.write(r"""\hline
\end{tabular}
""")
outfile.close()

## The following loop produces the 2 lines (4 entries) in the table
s = 4
for entry in [ # check the translation + flip groups
            (
            lambda m,n: ("1", m, n, s),
            lambda m,n: ("1", m, n, s+n),
            lambda m,n: FDR(1,1), # equal
            lambda m,n: n > m
        ),
            (
            lambda m,n: ("1", m, n, s),
            lambda m,n: ("1", m, n, -m-s),
                lambda m,n: "MIRROR",
            lambda m,n: n > m
        ),
            (
            lambda m,n: (".", m, n, s),
            lambda m,n: (".", m, n, s+n),
                lambda m,n: FDR(1,1), # equal
            lambda m,n: n > m
        ),
            (
            lambda m,n: (".", m, n, s),
            lambda m,n: (".", m, n, -m-s),
                lambda m,n: "MIRROR",
            lambda m,n: n > m
        ),]:
    print(f"check with parameter {s=}:",process_entry(entry))

print(f"OVERALL RUNTIME: {time() - start_overall}")
