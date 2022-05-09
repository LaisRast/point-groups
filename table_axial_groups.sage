#!/usr/bin/env python 
# coding: utf-8

## This script generates 4 LaTeX tables:
## * `tab_axial_combined.tex`:
##   For internal use, with the full information.
## * `tab_axial_full.tex`:
##   All axial groups with cross-references.
## * `tab_axial_pap.tex`:
##   Pyramidal and prismatic.
## * `tab_axial_hybrid.tex`:
##   Hybrid groups with methods.
## This files loads `table_polyhedral_groups.sage` to generate a catalog.
##
## Running time (with `precomputed_catalog`): 609.758593082428s

## Pass the following variables to `table_polyhedral_groups.sage`
only_axial = True
no_latex_file = True
load("table_polyhedral_groups.sage")

## Timing
start_time_axial = time()

## Functions
def is_axial(g):
    return all(e.act(Q(1)) in [Q(1), Q(-1)] for e in g)

def G3_projection(g):
    return set(
        e if e.act(Q(1))==Q(1)
        else e*FDR(Q(1),Q(1),True)
        for e in g)

def G3_plus_x4(g):
    return [e for e in g if e.act(Q(1))==Q(1)]

def classify_axial(g):
    g3 = G3_projection(g)
    g3plus = G3_plus_x4(g)
    if len(g3)==len(g3plus):
        assert set(g3)==set(g3plus),"55"
        if len(g3)==len(g):
            return "pyramidal"
        assert len(g)==2*len(g3), "66"
        return "prismatic"
    assert len(g)==len(g3)==2*len(g3plus),"44"
    return "hybrid axial"

## output dir
output_dir = "tables/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

######################
## Axial Full Table ##
######################
axial_combined_file = open(output_dir + "tab_axial_combined.tex","w")
axial_combined_file.write(r"""
\documentclass{article}
\usepackage[a4paper,margin=1cm]{geometry}
\usepackage{array}
\setlength{\extrarowheight}{1.8pt}
\begin{document}
\begin{table}
    \setlength{\extrarowheight}{1.8pt}
  \centering
%  \null\hskip-2,5cm
  \input tab_axial_pap
\end{table}
\begin{table}
    \setlength{\extrarowheight}{1.8pt}
  \centering
%  \null\hskip-2,5cm
  \input tab_axial_full
\end{table}

\begin{table}
    \setlength{\extrarowheight}{1.8pt}
    \centering
    % \null\hskip-2,5cm
    \input tab_axial_hybrid.tex
\end{table}

\end{document}
""")
axial_combined_file.close()
print("\nfile tab_axial_combined.tex written")

#########################################
## Axial Pyramidal and Prismatic Table ##
#########################################
axial_pap_file = open(output_dir + "tab_axial_pap.tex","w")
axial_pap_file.write(r"""
  \begin{tabular}[t]{|ccc|llr|llr|}
    \hline
    \multicolumn 3{|c|}{$G_3$}
&   \multicolumn 3{c|}{pyramidal groups $G_3 \times \{+x_4\}$}
&   \multicolumn 3{c|}{prismatic groups  $G_3 \times \{+x_4,-x_4\}$}
\\    \hline
    name % \cite{CS}
    & \hbox to 0pt{\hss orbitope\hss}& I.T.& CS name & Cox.%eter name
          & \llap{ord}er & CS  name & Cox.%eter
                                      \ name
          & \llap{o}rder \\\hline
""")
ii = 0
for g3 in poly3:
    pname,orbifold,IT,order = poly3groups[g3]
    axial_pap_file.write(fr"${pname}$& $\mathbf{{{orbifold}}}$& ${IT}$%"+"\n")

    pyrCS,pyrCox = pyramidal[g3]
    prisCS,prisCox = prismatic[g3]
    pyrname = CS_name_to_latex(pyrCS)
    prisname = CS_name_to_latex(prisCS)

    axial_pap_file.write(fr"&${pyrname}$& ${pyrCox}$ & ${order}$" )
    axial_pap_file.write(fr"&${prisname}$& ${prisCox}$ & ${order*2}$\\")
    if ii == len(poly3) - 1:
        axial_pap_file.write("[2pt]\n")
    else:
        axial_pap_file.write("\n")
    ii += 1 

axial_pap_file.write(r"""\hline
\end{tabular}
""")
axial_pap_file.close()
print("\nfile tab_axial_pap.tex written")

######################
## Axial Full Table ##
######################
def format_line(CSname, Coxname, desired_order, last_line):
    DV_number,DV_name,gens,order2,hurley = cat_axial[CSname]
    assert desired_order==order2, (desired_order,order2,CSname)
    CSlatex = CS_name_to_latex(CSname)
    if gens[0]==";":
        gens = gens[1:]
    g2 = gens_to_latex(gens)
    d = DV_number_format(DV_number)      
    output_line = (fr"&${CSlatex}$& {d}. ${DV_name}$&  ${Coxname}$ & {hurley}"+
                   fr"%& ${g2}$" + "\n" + fr"& ${order2}$\\")   
    if last_line:
        output_line += "[2pt]\n"
    else: 
        output_line += "\n"
    return output_line

cat_axial = {}
for (DV_number, DV_name, gens_CS, gens_RR), (order, name, hurley, diff) in zip(todo,out_list):
    if hurley is None:
        hurley = "n.cryst."
    else:
        hurley,page = BBNWZ_cat[hurley]
    cat_axial[name] = DV_number, DV_name, gens_RR, order, hurley

axial_full_file = open(output_dir + "tab_axial_full.tex","w")
header = r"""CS name & Du Val \# and name& Cox.%eter name 
 & {BBNZW} %& generators 
 & \llap{order} \\
\hline
"""
axial_full_file.write(r"""  \begin{tabular}[t]{|l|llllr|}
    \hline
 \multicolumn 6{|c|}{The 21 axial groups}
\\\hline
 \multicolumn 6{|c|}{pyramidal groups $G_3 \times \{+x_4\}$}
\\\hline $G_3$  &"""+header)
for ii in range(len(poly3)):
    g3  = poly3[ii]
    last_line = True if ii == len(poly3)-1 else False
    pname,orbifold,IT,order = poly3groups[g3]
    CS,Cox = pyramidal[g3]
    axial_full_file.write(fr"${pname}$")
    axial_full_file.write(format_line(CS, Cox, order, last_line))
axial_full_file.write(r"""\hline
 \multicolumn 6{|c|}{prismatic groups  $G_3 \times \{+x_4,-x_4\}$}
\\\hline $G_3$  &""" + header)
for ii in range(len(poly3)):
    g3  = poly3[ii]
    last_line = True if ii == len(poly3)-1 else False
    pname,orbifold,IT,order = poly3groups[g3]

    CS,Cox = prismatic[g3]
    axial_full_file.write(fr"${pname}$")
    axial_full_file.write(format_line(CS,Cox,2*order, last_line))
axial_full_file.write(r"""\hline
 \multicolumn 6{|c|}{hybrid axial groups $G_3^{+x_4}$ in $G_3$}
\\\hline $G_3^{+x_4}$ in $G_3$  &"""+header)
ii = 0
last_line = False
for (g3a,g3b),(CS,Cox,method) in hybrid.items():
    pnamea,_,_,ordera = poly3groups[g3a]
    pnameb,_,_,orderb = poly3groups[g3b]
    assert 2*ordera==orderb, (ordera,orderb,g3a,g3b)
    if ii == len(hybrid)-1: last_line = True
    axial_full_file.write(fr"${pnamea}$ in ${pnameb}$")
    axial_full_file.write(format_line(CS, Cox, 2*ordera, last_line))
    ii += 1

axial_full_file.write(r"""\hline
\end{tabular}
""")
axial_full_file.close()
print("\nfile tab_axial_full.tex written")

########################
## Axial Hybrid table ##
########################
axial_hybrid_file = open(output_dir + "tab_axial_hybrid.tex","w")
axial_hybrid_file.write(r"""\begin{tabular}[t]
%  {|l|lllr|l|}\hline
  {|l|llr|l|}\hline
    \multicolumn 5{|c|}{hybrid axial groups}\\
    \hline
    $G_3^{+x_4}$ in $G_3$ &    %$G(A;B)$
                            CS name %& Du Val \# and name
    & Coxeter name 
    & order & methods
    \\\hline""")

ii = 0
for (gr3_x4, gr3), (CS_name, Cox_name, methods) in hybrid.items():
    DuVal_number = DV_number_format(cat_axial[CS_name][0])
    DuVal_name = cat_axial[CS_name][1]
    order = cat_axial[CS_name][3]
    CS_name_latex = CS_name_to_latex(CS_name)

    gr3_x4_latex = poly3groups[gr3_x4][0]
    gr3_latex = poly3groups[gr3][0]

    assert order == poly3groups[gr3][3]

    axial_hybrid_file.write(f"${gr3_x4_latex}$ in ${gr3_latex}$ &${CS_name_latex}$ %& {DuVal_number} ${DuVal_name}$\n")
    axial_hybrid_file.write(f"&${Cox_name}$ & ${order}$ & {methods}\\\\")

    if ii == len(hybrid) - 1:
        axial_hybrid_file.write("[2pt]\n")
    else:
        axial_hybrid_file.write("\n")
    ii += 1
axial_hybrid_file.write("\\hline\n")
axial_hybrid_file.write("\\end{tabular}%")
axial_hybrid_file.close()

print(f"finished table_axial_groups.sage in: {time()-start_time_axial}s")
