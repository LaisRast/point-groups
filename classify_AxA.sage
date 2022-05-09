#!/usr/bin/env python 
# coding: utf-8

## Classifies all subgroups of $G=\pm[A\times A]$, 
## and their achiral extensions.
## The group $G$ is defined by the variable `name`.
## Output is written to the file `logs/{name}.log`.
## To direct output to `stdout`, change the variable `file` to `sys.stdout`.

load("main_functions.sage")

## input
name = "TxT" # or OxO or IxI
should_check_duplications = True

## functions
def achiral_extensions(AxA, name=""):
    ## Need cleaning
    """
    Return two catalogs:
    * The first one contains all subgroups of AxA group. These are chiral and
    presented as permutation groups.
    * The second one contains all possible achiral extensions of the subgroups
    in the first catalog. These are presented as AxB_group's
    """
    catalog_chiral = {}
    catalog_achiral = {}
    for n, g1 in enumerate(AxA.subgroups[:]):
        g1.AxB = AxA.convert_Pgroup_to_AxB(g1)
        clas = classify_group(g1.AxB)
        catalog_chiral.setdefault(clas,[]).append(g1)
        AxA.find_kernels(g1)
        classes = DisjointSet(AxA.Q1.q)
        for cose in cosets(AxA.Q1.q, g1.L0):
            rep = cose[0]
            for other in cose[1:]:
                classes.add(rep, other)

        reflections = set(classes.leader[c] for c in AxA.find_equal_pairs(g1))
        reflnew = []
        for ref in reflections: # preferred representative
                if AxA.Q1.one in classes.cluster[ref]:
                    ref = AxA.Q1.one
                elif AxA.Q1.minus_1 in classes.cluster[ref]:
                    ref = AxA.Q1.minus_1
                reflnew.append(ref)
        reflections = reflnew

        gens2 = [AxA.element_as_pair(e) for e in g1.gens()]
        for c in reflections:
            refl = FDR(Q(1),c,True)
            # check if normalizes
            g2 = set(refl.inverse() * e * refl for e in gens2)
            if g2 <= set(g1.AxB):
                newg = AxB_group(gens2+[refl])
                clas = classify_group(newg.group)
                catalog_achiral.setdefault(clas,[]).append((c,newg))

    return catalog_chiral, catalog_achiral

def check_duplicates(catalog, verbose=False, report_conjugacies=False, file=sys.stdout):
    """
    Check if groups in each category are
    a) not identical (should be by construction, in this case)
    b) isomorphic (this is a very fast test)
    c) even conjugate (dumb and slow)

    """  
    def get_AxB_group(gr):
        """
        gr can be a permutation group, a tuple (c, AxB_group) or an AxB_group.
        """
        if isinstance(gr, PermutationGroup_subgroup):
            gens = group = list(gr.AxB)
            new_gr = AxB_group(gens=gens, group=group)
            new_gr.original_pgroup = gr
            return new_gr
        elif isinstance(gr, tuple):
            assert isinstance(gr[1], AxB_group)
            return gr[1]
        elif isinstance(gr, AxB_group):
            return gr
        else:
            raise ValueError

    global DID_NOT_FIND
    DID_NOT_FIND = []

    nchecked = 0
    for cat, gs in catalog.items():
        if len(gs) > 1:

            ## convert to AxB group
            gs = [get_AxB_group(g) for g in gs]
            if verbose:
                print(f"# checking group fingerprint: {cat}", flush=True, file=file)
                print(f"order: {gs[0].order}", flush=True, file=file)
                print(f"number of groups: {len(gs)}", flush=True, file=file)

            ## check for being the identical
            if len(set(group_string_rep(g) for g in gs)) < len(gs):
                print("Warning: *identical* groups within", cat, flush=True, file=file)
                if None not in [g.name for g in gs]:
                    print("Among the following groups:", [g.name for g in gs], flush=True, file=file)
                # raise
            if verbose: print("identical test done", flush=True, file=file)

            ## check for isomophism
            nchecked_iso = 0
            for g2_index in range(1, len(gs)):
                g2 = gs[g2_index]
                g2_name = g2.name if g2.name is not None else f"gr_in_cat_index{g2_index}"
                g0 = gs[0]
                g0_name = g0.name if g0.name is not None else f"gr_in_cat_index{0}"
                nchecked_iso +=1
                if hasattr(g0, "original_pgroup") and hasattr(g2, "original_pgroup"):
                    isom = g0.original_pgroup.is_isomorphic(g2.original_pgroup)
                else:
                    isom = g2.is_isomorphic(g0)
                if not isom:
                    print("Warning: *nonisomorphic* groups within", cat, flush=True, file=file)
                    print(g0_name, "is not isomorphic to", g2_name, flush=True, file=file)
                    # raise
            if verbose: print(f"isomophism test done. {nchecked_iso} checks were performed", flush=True, file=file)

            ## check for conjugation
            nchecked_conj = 0
            for g2_index in range(1, len(gs)):
                g2 = gs[g2_index]
                g2_name = g2.name if g2.name is not None else f"gr_in_cat_index{g2_index}"
                g0 = gs[0]
                g0_name = g0.name if g0.name is not None else f"gr_in_cat_index{0}"
                nchecked_conj +=1
                try:
                    if g2.find_conjugation(g0): 
                        if report_conjugacies:
                            print(g2_name, "->", g0_name, ": conjugation by", g2.conj_name, flush=True, file=file)
                    else:
                        print("Warning: *nonconjugate* groups within", cat, flush=True, file=file)
                        print(g0_name, "is not conjugate to", g2_name, flush=True, file=file)
                except NotImplementedError:
                    print("Warning: *no conjugation found for* groups within", cat, flush=True, file=file)
                    print("The groups are", g0_name, "and", g2_name, flush=True, file=file)
                    DID_NOT_FIND.append([cat, g0, g2])
                    # raise
            if verbose: print(f"conjugacy test done. {nchecked_conj} checks were performed", flush=True, file=file)
            nchecked += len(gs)
    if verbose:
        print(f"duplication check finished. {nchecked} checks were performed.", flush=True, file=file)
    return DID_NOT_FIND

## set Q1 and Q2
if name == "TxT":
    Q1 = Q2 = group_2T
if name == "OxO":
    Q1 = Q2 = group_2O
if name == "IxI":
    Q1 = Q2 = group_2I

## output dir
Path("./logs/").mkdir(parents=True, exist_ok=True)

## output file
file = open(f"./logs/{name}.log", "w") # or sys.stdout

## run classification
start_time = time()
print(f"== classifying {name} ==", flush=True, file=file)

## generate the group
print("\n## generating group", flush=True, file=file)
time0()
gr = time1(AxB_product_group(Q1, Q2, name=name), file=file);
print("group order:", gr.order/2, flush=True, file=file)

## compute subgroups
print("\n## computing subgroups", flush=True, file=file)
time0()
time1(gr.compute_subgroups(only_full_groups=True), file=file);
print(len(gr.sgr0),"subgroups in total", flush=True, file=file)
print(len(gr.sgr1),"subgroups of them contain (-1,-1)", flush=True, file=file)

## print subgroups
# gr.print_subgroup_lengths()
lens = sorted({g.order()/2 for g in gr.subgroups})
print("subgroup orders:", lens, flush=True, file=file)
print("multiplicities: ", list((l,len(list(x for x in gr.subgroups if x.order()/2==l))) for l in lens), flush=True, file=file)

## compute achiral extensions
print("\n## computing achiral extensions", flush=True, file=file)
time0()
gr_cat = time1(achiral_extensions(gr,"gr"), file=file)

## print chiral groups
print("#catalog chiral:", flush=True, file=file)
for cl, grs in gr_cat[0].items():
    gr1 = grs[0]
    print(name, f"{len(grs):2}x  order", "%3s"%(gr1.order()/2), flush=True, file=file)

## print their achiral groups (select nice extension elements)
print("\n#catalog achiral:", flush=True, file=file)
for cl, grs in gr_cat[1].items(): 
    for c, gr1 in grs: 
        if str(c) == "1": 
            break 
    else: 
        for c, gr1 in grs: 
            if str(c)=="-1": 
                break 
        else: 
            c, gr1 = grs[0] 
    print (f"{name} {len(grs):2}x  order", "%3s"%(gr1.order/2), f'*[1,{c}]', flush=True, file=file)

if should_check_duplications:
    ## check duplications chiral
    print("\n## check duplications chiral", flush=True, file=file)
    time0()
    no_conjugation_found1 = time1(check_duplicates(gr_cat[0], verbose=True, report_conjugacies=True, file=file), file=file)

    if len(no_conjugation_found1):
        print("\n## did not find conjugation for the following chiral groups", flush=True, file=file)
        for cat, g0, g2 in no_conjugation_found1:
            print(cat, flush=True, file=file)
            print(g0, flush=True, file=file)
            print(g2, flush=True, file=file)
            print("", flush=True, file=file)
    else:
        print("\n## conjugations found for all the chiral groups", flush=True, file=file)


    ## check duplications achiral
    print()
    print("\n## check duplications achiral", flush=True, file=file)
    time0()
    no_conjugation_found2 = time1(check_duplicates(gr_cat[1], verbose=True, report_conjugacies=True, file=file), file=file)

    if len(no_conjugation_found2):
        print("\n## did not find conjugation for the following achiral groups", flush=True, file=file)
        for cat, g0, g2 in no_conjugation_found2:
            print(cat, flush=True, file=file)
            print(g0, flush=True, file=file)
            print(g2, flush=True, file=file)
            print("", flush=True, file=file)
    else:
        print("\n## conjugations found for all the achiral groups", flush=True, file=file)

print(f"\n## TOTAL TIME {time()-start_time:f}s", flush=True, file=file)
