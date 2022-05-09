#!/usr/bin/env python 
# coding: utf-8

## Generates order catalogs and stores them in `precomputed/` as json files.
## The keys are orders of the groups.
## The values are the names of the groups.
## Available counting functions are:
## * `count_all`
## * `count_tubular`
## * `count_chiral_toroidal`
## * `count_achiral_toroidal`
## * `count_polyhedral_and_axial`

load("main_functions.sage")

def fip_order(fip):
    """
    Return the order of the point group from an AxB group fingerprint.
    """
    order = sum([int(part.split(":")[1]) for part in fip.split()])
    assert is_even(order)
    return order/2
    
def toroidal_group_name(category, p1=None, p2=None, p3=None):
    """
    Return order, name and latex name of a toroidal group.
    """
    assert category in "1./\\X|+L*", "unknown category"

    if category == "1":
        m, n, s = p1, p2, p3
        name=f"G1^({s})_{m},{n}"
        latex_placeholder = r'\grp'+category+'_{{{1},{2}}}^{{({3})}}'
        if (n == 1) or (n == 2 and m % 2 == 1):
            latex_placeholder = r'\grp'+category+'_{{{1},{2}}}'
        order = m*n

    elif category == ".":
        m, n, s = p1, p2, p3
        name=f"G.^({s})_{m},{n}"
        latex_placeholder = r'\grp'+category+'_{{{1},{2}}}^{{({3})}}'
        if (n == 1) or (n == 2 and m % 2 == 1):
            latex_placeholder = r'\grp'+category+'_{{{1},{2}}}'
        order = 2*m*n

    elif category in "\\/":
        wallpaper, m, n = p1, p2, p3
        if wallpaper == "pm":
            assert m % 2 == n % 2 == 0
            multiplier = 1
        elif wallpaper == "pg":
            assert m % 2 == n % 2 == 0
            multiplier = 1
        elif wallpaper == "cm":
            assert (m - n) % 2 == 0
            multiplier = 2
        else:
            raise ValueError(wallpaper)
        name=f"G{category}^{wallpaper}_{m},{n}"
        if category == "/":
            latex_placeholder = r'\grp/'
        else:
            latex_placeholder = r'\grp\setminus'
        latex_placeholder += r'^{{\mathbf{{'+wallpaper+r'}}}}_{{{2},{3}}}'
        order = multiplier*m*n

    elif category == "X":
        wallpaper, m, n = p1, p2, p3
        if wallpaper == "cmm":
            assert (m - n) % 2 == 0
            multiplier = 4
        else:
            assert m % 2 == n % 2 == 0
            multiplier = 2
        name=f"G{category}^{wallpaper}_{m},{n}"
        wallpaper_full = wallpaper[:1]+"2"+wallpaper[1:]  # according to new IT
        latex_placeholder = r'\grp X^{{\mathbf{{'+wallpaper_full+r'}}}}_{{{2},{3}}}'
        order = multiplier*m*n

    elif category == "|":
        wallpaper, m, n = p1, p2, p3
        if wallpaper == "pm":
            multiplier = 2
        elif wallpaper == "cm":
            multiplier = 4
        elif wallpaper == "pg":
            multiplier = 2
        else:
            raise ValueError(wallpaper)
        name=f"G{category}^{wallpaper}_{m},{n}"
        latex_placeholder = r'\grp'+category+ r'^{{\mathbf{{'+wallpaper+r'}}}}_{{{2},{3}}}'
        order = multiplier*m*n
    
    elif category == "+":
        wallpaper, m, n = p1, p2, p3
        if wallpaper == "pmm":
            multiplier = 4
        elif wallpaper == "pmg":    # mirror is vertical
            multiplier = 4
        elif wallpaper == "pgg":
            multiplier = 4
        elif wallpaper == "cmm":
            multiplier = 8
        name=f"G{category}^{wallpaper}_{m},{n}"
        wallpaper_full = wallpaper[:1]+"2"+wallpaper[1:] # according to new IT
        latex_placeholder = r'\grp'+category+ r'^{{\mathbf{{'+wallpaper_full+r'}}}}_{{{2},{3}}}'
        order = multiplier*m*n

    elif category == "L":
        a, b = p1, p2
        assert p3 is None
        wallpaper = "p4"
        name=f"G{category}_{a},{b}"
        order = 4*(a^2+b^2)
        latex_placeholder = r'\grp L_{{{1},{2}}}'

    elif category == "*":
        wallpaper, n = p1, p2
        assert p3 is None
        if wallpaper in ["p4mS", "p4gS"]:
            multiplier = 16
        elif wallpaper in ["p4mU", "p4gU"]:
            multiplier = 8
        else:
            raise ValueError(wallpaper)
        name=f"G{category}^{wallpaper[2:]}_{n}" # cut out "p4"
        order = multiplier*n*n
        wallpaper_full = (r'\mathbf{{'+wallpaper[:3]
                          +'m'     # according to new IT
                          +r"}}\textrm{{"+wallpaper[-1]+"}}")
        latex_placeholder = r'\grp'+category+ r'^{{'+wallpaper_full+r'}}_{{{2}}}'

    latex = latex_placeholder.format(0,p1,p2,p3)
    return order, name, latex

def count_chiral_toroidal(orderbound, latex_only=False):
    """
    Return the number of the chiral toroidal groups up to `orderbound` and a
    catalog of them grouped by order.
    """

    order_catalog = {}
    num_groups = 0

    ## grp. and grp1
    for typ in "1.":
        for m in range(1, 1+orderbound):
            for n in range(1, 1+orderbound):
                if m*n > orderbound or (typ=="." and 2*m*n>orderbound): break
                for s in range(-(m//2),1+(n-m)//2):
                    if typ=="." and m<=2 and n==1:
                        continue 
                    order, name, latex = toroidal_group_name(typ,m,n,s)
                    if latex_only: name = latex
                    order_catalog.setdefault(order, []).append(name)
                    num_groups += 1

    ## grp\ and grp/
    for m in range(1,1+orderbound//2):
        for n in range(1,1+orderbound//2):
            if 2*m*n>orderbound: break
            for typ in "\\/":
                for wallp in ("pm","pg","cm"):
                    if wallp != "cm":
                        if 4*m*n > orderbound: continue
                        order, name, latex = toroidal_group_name(typ,wallp,2*m,2*n)
                    if wallp == "cm":
                        if (m-n)%2 != 0: continue
                        order, name, latex = toroidal_group_name(typ,wallp,m,n)
                    if ((wallp in ("pm","cm") and (n==1 or m==1))
                        or (wallp=="pg" and ((typ=="/" and n==1) or (typ=="\\" and m==1)))
                        or (wallp=="cm" and ((typ=="/" and n==2) or (typ=="\\" and m==2)))
                        ): 
                        continue 
                    if latex_only: name = latex
                    order_catalog.setdefault(order, []).append(name)
                    num_groups += 1

    ## grpX
    for m in range(1,1+orderbound//4):
        for n in range(1,1+orderbound//4):
            if 4*m*n>orderbound: break
            for wallp in ("pmm","pgg","cmm","pmg","pgm"):
                    if wallp == "cmm":
                        if (m-n)%2 != 0: continue
                        order, name, latex = toroidal_group_name("X",wallp,m,n)
                    else:
                        if 8*m*n>orderbound: continue
                        order, name, latex = toroidal_group_name("X",wallp,2*m,2*n)
                    if n==1 or m==1 or (wallp=="cmm" and (n==2 or m==2)):                     
                        continue 
                    if latex_only: name = latex
                    order_catalog.setdefault(order, []).append(name)
                    num_groups += 1
    
    assert len(set(flatten((list(order_catalog.values()))))) == num_groups
    return num_groups, order_catalog

def count_achiral_toroidal(orderbound, latex_only=False):
    """
    Return the number of the achiral toroidal groups up to `orderbound` and a
    catalog of them grouped by order.
    """

    order_catalog = {}
    num_groups = 0

    ## grp|
    for wallp in ("pm", "cm", "pg"):
        for m in IntegerRange(1, orderbound+1):
            for n in IntegerRange(1,orderbound+1):
                if wallp == "cm" and 4*m*n > orderbound: break
                else:
                    if 2*m*n > orderbound: break
                order, name, latex = toroidal_group_name("|", wallp, m, n)
                if latex_only: name = latex
                order_catalog.setdefault(order, []).append(name)
                num_groups += 1
   
    ## grp+
    for wallp in ("pmm", "pmg", "pgg", "cmm"):
        for m in IntegerRange(1, orderbound+1):
            for n in IntegerRange(1, orderbound+1):
                if wallp == "cmm" and 8*m*n > orderbound: break
                else:
                    if 4*m*n > orderbound: break
                if (m == 1 and n == 1) or ((m < n) and (wallp !="pmg")):
                    continue
                order, name, latex = toroidal_group_name("+", wallp, m, n)
                if latex_only: name = latex
                order_catalog.setdefault(order, []).append(name)
                num_groups += 1
   
    ## grpL
    for a in IntegerRange(1,orderbound+1):
        for b in IntegerRange(0,orderbound+1):
            if 4*(a^2+b^2) > orderbound: break
            if ((a,b) in [(1,0), (1,1), (2,0)]) or (a < b):
                continue
            order, name, latex = toroidal_group_name("L", a, b)
            if latex_only: name = latex
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1
   
    ## grp*
    for wallp in ("p4mU", "p4gU", "p4mS", "p4gS"):
        for n in IntegerRange(1,orderbound+1):
            if wallp in ("p4mS", "p4gS") and 16*n*n > orderbound: break
            else:
                if 8*n*n > orderbound: break
            if (wallp in ("p4mS", "p4gS") and n == 1) or (wallp in ("p4mU", "p4gU") and (n == 1 or n == 2)):
                continue
            order, name, latex = toroidal_group_name("*", wallp, n)
            if latex_only: name = latex
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    assert len(set(flatten((list(order_catalog.values()))))) == num_groups
    return num_groups, order_catalog

def test_count_toroidal():
    """
    Check outputs of `count_chiral_toroidal` and `count_achiral_toroidal`
    using precomputed data.
    """
    ## chiral toroidal
    with open("./precomputed/toroidal_chiral_catalog_orderbound_130.json") as f:
        catalog_precomputed = json.load(f)
    order_catalog = count_chiral_toroidal(130)
    for fip, grs in catalog_precomputed.items():
        order = fip_order(fip)
        assert order in order_catalog, f"order_catalog is missing groups of order {order}: {' '.join(grs)}"
        for gr in grs:
            assert gr in order_catalog[order], f"order_catalog is missing group of order {order}: {gr}"
    
    order_catalog_set = set(flatten((list(order_catalog.values()))))
    catalog_precomuted_set = set(flatten((list(catalog_precomputed.values()))))
    assert order_catalog_set == catalog_precomuted_set, "order_catalog has groups not in precomputed_catalog"

    ## achiral toroidal
    with open("./precomputed/toroidal_achiral_catalog_orderbound_130.json") as f:
        catalog_precomputed = json.load(f)
    order_catalog = count_achiral_toroidal(130)
    for fip, grs in catalog_precomputed.items():
        order = fip_order(fip)
        assert order in order_catalog, f"order_catalog is missing groups of order {order}: {' '.join(grs)}"
        for gr in grs:
            assert gr in order_catalog[order], f"order_catalog is missing group of order {order}: {gr}"
    
    order_catalog_set = set(flatten((list(order_catalog.values()))))
    catalog_precomuted_set = set(flatten((list(catalog_precomputed.values()))))
    assert order_catalog_set == catalog_precomuted_set, "order_catalog has groups not in precomputed_catalog"
    return

def count_polyhedral_and_axial(typ, orderbound, latex_only=False):
    """
    Return the number of the polyhedral and the axial groups up to
    `orderbound` and a catalog of them grouped by order.
    Here `typ` can be `chiral` or `achiral` or `all`.
    """

    order_catalog = {}
    num_groups = 0

    with open("./precomputed/polyhedral_and_axial_catalog.json") as f:
        catalog_precomputed = json.load(f)
    for fip, grs in catalog_precomputed.items():
        if typ == "chiral" and "*" in fip: continue
        if typ == "achiral" and "*" not in fip: continue
        order = fip_order(fip)
        if order > orderbound: continue
        for name in grs:
            if latex_only: name = CS_name_to_latex(name)
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    assert len(set(flatten((list(order_catalog.values()))))) == num_groups
    return num_groups, order_catalog


def count_tubular(orderbound, latex_only=False):
    """
    Return the number of the tubular groups up to `orderbound` and a catalog
    of them grouped by order.
    """

    order_catalog = {}
    num_groups = 0

    for (gr3_name, gr3_order) in [("I", 120), ("O", 48), ("T", 24)]:
        ## +-[LxCn]
        for n in range(1, orderbound//gr3_order+1): # 2 because of +-
            for name1, name2 in [(f"C{n}", gr3_name),
                                 (gr3_name, f"C{n}")]:
                name = f"+-[{name1}x{name2}]"
                if latex_only: name = CS_name_to_latex(name)
                order = gr3_order*n
                order_catalog.setdefault(order, []).append(name)
                num_groups += 1

        ## +-[LxD2n]
        for n in range(2, orderbound//(2*gr3_order)+1):
            for name1, name2 in [(f"D{2*n}", gr3_name),
                                 (gr3_name, f"D{2*n}")]:
                name = f"+-[{name1}x{name2}]"
                if latex_only: name = CS_name_to_latex(name)
                order = 2*gr3_order*n
                order_catalog.setdefault(order, []).append(name)
                num_groups += 1

    ## +-1/2[OxD2n]
    factor = 48
    for n in range(2, orderbound//factor+1): #n=1 : D2 becomes C2
        for name1, name2 in [("O", f"D{2*n}"),
                             (f"D{2*n}", "O")]:
            name = f"+-1/2[{name1}x{name2}]"
            if latex_only: name = CS_name_to_latex(name)
            order = 48*n
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    ## +-1/2[OxD'4n]
    factor = 96        
    for n in range(2, orderbound//factor+1): #n=1: Dbar4 becomes D4
        for name1, name2 in [("O", f"D'{4*n}"),
                             (f"D'{4*n}", "O")]:
            name = f"+-1/2[{name1}x{name2}]"
            if latex_only: name = CS_name_to_latex(name)
            order = 96*n
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    ## +-1/6[OxD6n]
    factor = 48
    for n in range(1, orderbound//factor+1):
        for name1, name2 in [("O", f"D{6*n}"),
                             (f"D{6*n}", "O")]:
            name = f"+-1/6[{name1}x{name2}]"
            if latex_only: name = CS_name_to_latex(name)
            order = 48*n
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1
 
    ## +-1/2[OxC2n]
    factor = 48
    for n in range(1, orderbound//factor+1):
        for name1, name2 in [("O", f"C{2*n}"),
                             (f"C{2*n}", "O")]:
            name = f"+-1/2[{name1}x{name2}]"
            if latex_only: name = CS_name_to_latex(name)
            order = 48*n
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    ## +-1/3[TxC3n]
    factor = 24
    for n in range(1, orderbound//factor+1):
        for name1, name2 in [("T", f"C{3*n}"),
                             (f"C{3*n}", "T")]:
            name = f"+-1/3[{name1}x{name2}]"
            if latex_only: name = CS_name_to_latex(name)
            order = 24*n
            order_catalog.setdefault(order, []).append(name)
            num_groups += 1

    assert len(set(flatten((list(order_catalog.values()))))) == num_groups
    return num_groups, order_catalog

def test_count_tubular():
    """
    Check outputs of `count_tubular` using precomputed data.
    """
    with open("./precomputed/tubular_catalog_orderbound_2000.json") as f:
        catalog_precomputed = json.load(f)
    order_catalog = count_tubular(2000)
    for fip, grs in catalog_precomputed.items():
        order = fip_order(fip)
        assert order in order_catalog, f"order_catalog is missing groups of order {order}: {' '.join(grs)}"
        for gr in grs:
            assert gr in order_catalog[order], f"order_catalog is missing group of order {order}: {gr}"
    
    order_catalog_set = set(flatten((list(order_catalog.values()))))
    catalog_precomuted_set = set(flatten((list(catalog_precomputed.values()))))
    assert order_catalog_set == catalog_precomuted_set, "order_catalog has groups not in precomputed_catalog"
    return

def count_all(typ, orderbound, latex_only=False):
    """
    Return the number of point groups up to `orderbound` and a catalog of
    them grouped by order.
    Here `typ` can be `chiral` or `achiral` or `all`.
    """
    order_catalog = {}
    num_groups = 0
    
    order_catalogs = [count_polyhedral_and_axial(typ, orderbound, latex_only)]
    if typ == "chiral":
        order_catalogs.append(count_tubular(orderbound, latex_only))
        order_catalogs.append(count_chiral_toroidal(orderbound, latex_only))
    elif typ == "achiral":
        order_catalogs.append(count_achiral_toroidal(orderbound, latex_only))
    else:
        order_catalogs.append(count_tubular(orderbound, latex_only))
        order_catalogs.append(count_chiral_toroidal(orderbound, latex_only))
        order_catalogs.append(count_achiral_toroidal(orderbound, latex_only))

    for _, order_catalog2 in order_catalogs:
        for ordr, grs in order_catalog2.items():
            for gr in grs:
                order_catalog.setdefault(int(ordr), []).append(gr)
                num_groups += 1

    return num_groups, order_catalog

def main():
    ## input
    should_generate_count_all = True
    orderbound = 130
    for typ in ["all", "chiral", "achiral"]:

        ## output directory
        output_dir = "precomputed/"
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if should_generate_count_all:
            order_catalog = count_all(typ, orderbound, latex_only=False)[1]
            filename = f"order_catalog_{typ}_orderbound_{orderbound}.json"
            with open(output_dir + filename, "w") as f:
                json.dump(order_catalog, f, indent=4)

if __name__ == "__main__":
    main()
