#!/usr/bin/env python 
# coding: utf-8

## Classifies all chiral toroidal groups up to order `orderbound//2``,
## and the achiral toroidal groups up to order `orderbound`.
## Output is written to the file `logs/DxD_{orderbound}.log`.
## To direct output to `stdout`, change the variable `file` to `sys.stdout`.

load("main_functions.sage")
load("conjugacy_Dn.sage")

start_time = time()

## input
orderbound = 400
orderbound_table = orderbound

## output dir
Path("./logs/").mkdir(parents=True, exist_ok=True)

## output file
file = open(f"./logs/DxD_{orderbound}.log", "w") # or sys.stdout

## sage version 
print(sage.misc.banner.version(), file=file, flush=True)

print(f"elapsed time for loading: {time()-start_time}s", file=file, flush=True)

## increase GAP pre-set memory limit
if sage.misc.banner.require_version(major=9, minor=3):
  # sage.interfaces.gap.gap_cmd = 'gap -r -o 24G '
  sage.interfaces.gap.gap_cmd = 'gap -r -o 50G '
else:
    # The following works in sage 9.2, but no longer in sage 9.5:
    from sage.interfaces.gap import set_gap_memory_pool_size, get_gap_memory_pool_size
    print(f"GAP default memory pool size {get_gap_memory_pool_size()}", file=file, flush=True)
    set_gap_memory_pool_size(50* 10**9)
    print(f"GAP adjusted memory pool size {get_gap_memory_pool_size()}", file=file, flush=True)

## ensure enough memory for GAP
if 0:
    start_time192 = time()
    n = 192
    grCD = Dn_group(n)
    CDnxCDn2 = AxB_product_group(grCD,achiral=True,limit=8*n*n)
    print(f"Generated. {CDnxCDn2.order=} {CDnxCDn2.name=}", file=file) 
    print(f"elapsed time: {time()-start_time192}s", file=file)
 

## functions
def lookup_in_catalog(gr0):
    """
    Check if `gr0` is in catalog.
    Here `gr0` is a list of elements.
    """
    str0 = group_string_rep(gr0)
    clas0 = classify_group(gr0)
    
    visited.add(clas0)

    grs = catalog[clas0]

    ## check for identical
    for gr in grs:
        if str0 == group_string_rep(gr):
            return ("identical", gr)

    ## check for conjugacy
    for gr in grs:
        try:
            OK, congL, congR = find_conjugation_2D(gr.group, gr0, name1=gr.name)
            if OK:
                return ("conjugate to", gr, congL, congR)
        except ValueError:
            ## switch to old method
            if gr.find_conjugation(gr0):
                return ("conjugate to", gr)
    print(f"No conjugation found for gr0 = {str0}", flush=True, file=file)
    raise ValueError

## run classification
print(f"== classifying DxD.2 == {orderbound=}", flush=True, file=file)
print(f"Elapsed time for preparation: {time()-start_time}s\n",
      flush=True, file=file)

## generate catalog
print(f"\n== catalog of toroidal groups ==", flush=True, file=file)
try:
    print(f"{len(catalog)} groups exist in catalog", flush=True, file=file)
except NameError:
    catalog = {}
    generate_chiral_toroidal(orderbound=orderbound_table//2, generate_full=True, verbose=False)
    n_chiral = len(catalog)
    print(f"{n_chiral} chiral groups up to order {orderbound_table//2}", flush=True, file=file)
    generate_achiral_toroidal(orderbound=orderbound_table, generate_full=True, verbose=False)
    print(f"{len(catalog)-n_chiral} achiral groups up to order {orderbound_table}", flush=True, file=file)

## generate groups systematically
visited = set()
total_it = 0

## generate the groups systematically
for n in range(orderbound//2, orderbound//4, -1): 
    grs = []
    if n % 2 == 0:
        grs.append((Dn_group(n), "D"))
    grs.append((Cn_group(n), "C"))

    for grCD, name in grs:

        print(f"\n=== classifying {name}{n}x{name}{n}.2 ===", flush=True, file=file)
        gr_start_time = time()

        ## generate the group
        print("## generating group", flush=True, file=file)
        time0()
        gr = time1(AxB_product_group(grCD, achiral=True, limit=8*n*n), file=file)
        # print(f"name: {name}{n}x{name}{n}", flush=True, file=file)
        print(f"order: {gr.order/2}", flush=True, file=file)
        assert gr.order == 2*n*n*2*2

        ## compute subgroups
        print("\n## computing subgroups", flush=True, file=file)
        time0()
        time1(gr.compute_subgroups(only_full_groups=False), file=file)
        print(len(gr.sgr0),"subgroups in total", flush=True, file=file)
        # print(len(gr.sgr1),"subgroups of them contain (-1,-1)", flush=True, file=file)

        ## print subgroups
        # gr.print_subgroup_lengths()
        lens = sorted({g.order()/2 for g in gr.subgroups})
        print("subgroup orders:", lens, flush=True, file=file)
        print("multiplicities: ", list((l,len(list(x for x in gr.subgroups if x.order()/2==l))) for l in lens), flush=True, file=file)

        ## checking subgroups
        print("\n## checking subgroups", flush=True, file=file)
        it=0
        for it2, subgroup in enumerate(gr.subgroups):
            if subgroup.order() > orderbound_table:
                # Since gr is achiral, this is the true order
                continue
            gr0 = gr.convert_Pgroup_to_AxB(subgroup)
            achiral = any(e.star for e in gr0)
            if not achiral and len(gr0) > orderbound_table:
                continue
            result = lookup_in_catalog(gr0)
            it += 1
            if it % 1000==0:
                found_text = result[0]
                g = result[1]
                print(f"subgroup #{it}; order {subgroup.order()}; {found_text} {g.name}", flush=True, file=file)
        print(f"\n{name}{n}x{name}{n} done", flush=True, file=file)
        print(f"{it} groups were checked",flush=True, file=file)
        print(f"Elapsed time: {time()-gr_start_time}s\n", flush=True, file=file)
        total_it += it

print("finished checking subgroups", flush=True, file=file)

## missed groups
print("== missed ==", flush=True, file=file)
missed = set(catalog.keys()) - visited;
print (f"{len(missed)} groups missed, {len(visited)} groups visited", flush=True, file=file)
for cat in list(missed)[:0]:
    gr = catalog[cat][0]
    print("missed", gr.name, gr.order, flush=True, file=file)
    
## summary
print("== summary ==", flush=True, file=file)
print(f"{total_it} groups were checked in total", flush=True, file=file)
print(f"{len(visited)} groups visited", flush=True, file=file)
print(f"{len(missed)} groups missed", flush=True, file=file)
print(f"\n## TOTAL TIME {time()-start_time:f}s", flush=True, file=file)
