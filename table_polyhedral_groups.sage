#!/usr/bin/env python 
# coding: utf-8

## This script has two parts:
## * Checks that our generators and CS generators for the polyhedral and the
##   axial groups generate geometrically equivalent (but not necessarily equal)
##   groups.
##   The list of generators is given in `polyhedral_and_axial_lists.sage`.
## * Generates 3 LaTeX tables:
##     * `tab_polyhedral_and_axial_check.tex`:
##       For internal use, with the full information.
##     * `tab_polyhedral_and_axial_appendix.tex`:
##       All polyhedral and axial groups, with generators.
##     * `tab_polyhedral.tex`:
##       An overview of the polyhedral groups.
## Running time (with `precomputed_catalog`): 1631.6972992420197s

load("main_functions.sage")
load("BBNWZ_cat.sage")

## input
precomputed_catalog = True
generate_latex = True

## Timing
start_time_polyhedral = time()

## Initialize polyhedral and axial catlog
if precomputed_catalog:
    with open("precomputed/polyhedral_and_axial_catalog.json") as f:
        catalog = json.load(f)
else:
    catalog = dict()
    generate_polyhedral_and_axial(names_only=True)

## Load polyhedral and axial groups data (names, gens)
load("polyhedral_and_axial_lists.sage")
todo = polyhedral_list_Cox + axial_list

## Variables setten from outside this scripts
try:
    if only_axial:
        print("Only axial groups considered")
        todo = axial_list
except:
    pass

try:
    if no_latex_file:
        generate_latex = False
except:
    pass

## Get the Coxeter names for the axial groups
load("axial_data.sage")
axial_Cox_names = {it[0]: it[1] for it in list(pyramidal.values())+list(prismatic.values())+list(hybrid.values())}

## String to quaternion and latex table
gens_table = {
  # code: (FDR,                  latex),
    "-":  (FDR(Q(1),Q(-1)),      "-1"),
    "-*": (FDR(Q(1),Q(-1),True), "-{*}"), 
    "*":  (FDR(Q(1),Q(1),True),  "{*}"),
    ";":  (FDR(Q(1),Q(1)),        ";"), # identity
    }

## Functions
def translate_gens(e):
    if type(e) is str:
        return gens_table[e][0]
    star = e[0]=="*"
    if star:
        e = e[1:]
    return FDR(e[0],e[1],star)

def gens_to_latex(gens):
    comma = False
    result = ""
    for e in gens:
        if e==";":
            result += r"\textbf{;}"
            comma = False
        else:
            if comma:
                result += ","
            if type(e) is str:
                result += gens_table[e][1]
            else:
                star = e[0]=="*"
                if star:
                    result += "{*}"
                    e = e[1:]
                result += f"[{quaternion_to_latex(e[0])},{quaternion_to_latex(e[1])}]"
            comma = True
    return result

def DV_number_format(d):
    end = ""
    d = str(d)
    while d.endswith("'"):
        end += "'"
        d = d[:-1]
    if end:
        d += r"\rlap{$"+end+"$}"
    return d

def remove_blanks(s):
    "remove blanks that were inserted to influence the sorting order"
    return "".join(c for c in s if c!=" ")
 
## Generate and check groups
## It is okay if CS group and our group are not identical.
print('Part 1: Generate and check the groups')
out_list = []
for item in todo:
    DV_number, DV_name, gens_CS, gens_RR = item[:4]
    gens = [translate_gens(e) for e in gens_CS]
    gr = mk_group(gens, one=FDR(Q(1), Q(1)), limit=28800)
    assert FDR(Q(-1), Q(-1)) in gr, (DV_number, DV_name, "not closed under -1")
    clas = classify_group(gr)
    try:
        assert len(catalog[clas]) == 1
        data = len(gr)/2, catalog[clas][0], Hurley_pattern(clas)
    except KeyError:
        print( len(gr)/2, gens_CS)
        print( len(gr)/2, gens, DV_number)
        raise
    print(DV_number, *data, flush=True)
    gens_RR = [translate_gens(e) for e in gens_RR]
    gr_RR = mk_group(gens_RR, one=FDR(Q(1),Q(1)), limit=33333)
    if classify_group(gr_RR) != classify_group(gr):
        data += ("!!!!!!!!!!!!!",)
        print("\n  *****************groups are not the same!:", DV_number, DV_name,"\n")
        raise
    elif set(gr_RR) != set(gr):
        print("\n  *****************groups are not identical:", DV_number, DV_name,"\n")
        data += (" DIFF",)
    else:
        data += ("",)
    out_list.append(data)

if generate_latex:
    print('Part 2: Make the LaTeX table')
    ## output dir
    output_dir = "tables/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    outname = output_dir + "tab_polyhedral_and_axial_check.tex"
    outfile = open(outname,"w")
    outfile.write(r"""
    \documentclass{article}
    \usepackage[a4paper,landscape,margin=7mm,top=3mm]{geometry}
    \usepackage{array}
    \setlength{\extrarowheight}{0pt}%{0,5pt}
    \begin{document}
    \small
    \footnotesize
     DIFF means: the groups are not identical.

    """)
    tablehead=r"""
    \begin{tabular}{|llllllrl|}
      D\rlap{u Val \# and name} &&
      CS name& Cox.\ name & generators [CS] & generators (ours) & order & BBNWZ\#\\\hline
    """
    tablefoot =r"""[2pt]\hline
    \end{tabular}
    """

    lines = []
    lines_appendix = []
    for item, (order, name, hurley, diff) in zip(todo, out_list):
        DV_number, DV_name, gens_CS, gens_RR = item[:4]
        try:
            Cox_name = item[4]
        except IndexError:
            try:
                Cox_name = axial_Cox_names[name]
            except KeyError:
                Cox_name = "-"
        if hurley is None:
            hurley = "not cryst."
        else:
            hurley, page = BBNWZ_cat[hurley]
            # hurley = "".join("$'$" if c=="'" else "" if c=="*" else
            #                  r"\," if c==" " else c   for c in hurley)
            # if page!=100:
            #     hurley=BB
        g1 = gens_to_latex(gens_CS)
        g2 = gens_to_latex(gens_RR)
        lname = CS_name_to_latex(name)
        d = DV_number_format(DV_number)
        text = (fr'{d}.&${DV_name}$&${lname}$& ${Cox_name}$&${g1}${diff}& ${g2}$&{order}&{hurley}\\')
        lines.append(text +'\n')
        text = fr'{d}. ${DV_name}$&${lname}$& ${g2}$& ${Cox_name}$&{order}&{hurley}\\'
        lines_appendix.append(text +'\n')
        print(d, end=" ")
    print()
    outfile.write(tablehead)
    outfile.write("".join(lines))
    outfile.write(tablefoot)
    outfile.write(tablehead)
    outfile.write("".join(sorted(lines)))
    outfile.write(tablefoot)
    outfile.write(r"""
    \end{document}
    """)
    outfile.close()
    print("\nfile",outname,"written")



    outname = output_dir + "tab_polyhedral_and_axial_appendix.tex"
    outfile = open(outname, "w")
    outfile.write(r"""
    %\documentclass{article}
    %\usepackage[a4paper,margin=7mm]{geometry}
    %\usepackage{array}
    %\setlength{\extrarowheight}{0,6pt}
    %\begin{document}
    \begin{tabular}{|l@{ }l@{ \ }ll@{ }r@{ \ }l|}
      \hline
      D\rlap{u Val \# \& name} &
      CS name & generators& Cox.\ name & \llap{ord}er & BBNWZ\\\hline
    """)
    prev = ""
    for line in sorted(lines_appendix):
        if prev.startswith("39."):
            outfile.write(r"[2pt]\hline"+"\n")
        if remove_blanks(line)[:50] == remove_blanks(prev)[:50]:
            parts1 = prev.split("&")
            parts2 = line.split("&")
            assert parts1[3:6]==parts2[3:6]
            assert parts1[5].endswith(r'\\'+'\n')
            order = parts1[4]
            #if len(order)<3: order = r'\phantom0'+order
            neworder = r'\smash{\lower 1,4ex\hbox{'+order+'}}'
            newBBNZ =  r'\smash{\lower 1,4ex\hbox{'+parts1[5][:-3]+'}}' # remove "\\"
            newCox =   r'\smash{\lower 1,4ex\hbox{\llap{$\biggr\}$ }'+parts1[3].strip()+'}}'
            name1 = parts1[0].split()
            name2 = parts2[0].split()
            #name1 = r'\phantom{'+name1[0]+"}"
            ## sometimes, the names contain extra blanks, for sorting order
            
            # FIRST THE DU VAL (parts2), THEN THE COXETER-SMITH (parts1)
            prev = "&".join(
                [" ".join(name2),"",parts2[2], newCox,neworder,
                 newBBNZ + r'\\'+'\n'])
            line = "&".join(
                 ["", parts1[1], parts1[2],"","",
                 r'\\'+'\n'])
        outfile.write(prev)
        prev = line
    outfile.write(prev)
    outfile.write(tablefoot)
    outfile.write(r"%\end{document}"+"\n")
    outfile.close()
    print("\nfile",outname,"written")

    lines = []
    prev = ""
    for ((DV_number,DV_name,gens_CS,gens_RR,Cox_name,method),
         (order,name,hurley,diff)) in zip(polyhedral_list_Cox,out_list):
        if hurley is None:
            hurley = "n.cryst."
        else:
            hurley,page = BBNWZ_cat[hurley]
        if gens_RR[0]==";":
            gens_RR = gens_RR[1:]
        g2 = gens_to_latex(gens_RR)
        lname = CS_name_to_latex(name)
        d = DV_number_format(DV_number)
        # text = fr'${lname}$& {d}. ${DV_name}$& ${Cox_name}$& ${g2}$&{order}&{hurley}&{method}\\'
        text = fr'${lname}$& {d}. ${DV_name}$& ${Cox_name}$&{{{order}}}&{method}\\'
        print(d,end=" ")
        if remove_blanks(text) != remove_blanks(prev):
            lines.append(text+'\n')
        prev = text	
    outname = output_dir+"tab_polyhedral.tex"
    outfile=open(outname, "w")
    headline = r"[2pt]\hline \multicolumn 5{|l|}{symmetries of %s}\\[1pt]\hline"+"\n"
    outfile.write(r"""%\documentclass{article}
    %\begin{document}
    \begin{tabular}[t]{|lllr|l|}\hline
    CS name & Du Val \# and
              name& Coxeter name& \llap{ord}er &method\\
    """
    + headline % r"the 120-cell $Q_{123}=\{5,3,3\}$ / the 600-cell $P_{600}=\{3,3,5\}$")
    outfile.write("".join(lines[:6]))
    outfile.write(headline % r"the 24-cell $P_T=\{3,4,3\}$ and its polar 24-cell $P_{T_1}$")
    outfile.write("".join(lines[6:11]))
    outfile.write(r"[2pt]\hline"+"\n") # the first five analogous to the CUBE in 3d
    outfile.write("".join(lines[11:15]))
    outfile.write(headline % (r"the hypercube $\{4,3,3\}$ /" +
                              r" the cross-polytope $\{3,3,4\}$"))
    outfile.write("".join(lines[15:-5]))
    outfile.write(headline % r"the simplex $\{3,3,3\}$ and its polar")
    outfile.write("".join(lines[-5:]))
    outfile.write(r"[2pt]\hline \end{tabular}"+"\n")
    outfile.write(r"%\end{document}"+"\n")
    outfile.close()
    print("\nfile", outname, "written")

print(f"finished table_polyhedral_groups.sage in: {time()-start_time_polyhedral}s")
