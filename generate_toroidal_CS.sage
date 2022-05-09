#!/usr/bin/env python 
# coding: utf-8

## Contains functions to generate catalogs for the
## toroidal groups using CS classification:
## Available functions:
## * `generate_chiral_toroidal_CS`
## * `generate_achiral_toroidal_CS`
## * `make_toroidal_cross_reference`

def generate_chiral_toroidal_CS(orderbound=130, generate_full=False,
                                names_only=False, latex_only=False,
                                verbose=True):
    """
    Generate chiral toroidal groups (CS) and insert them in `catalog`.

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
    130        | True          | 28027      | 172.60427021980286s
    130        | False         | 9666       | 73.40611362457275s
    500        | True          | 407382     | 10267.683957099915s
    500        | False         | 134099     | 3921.1668684482574s
    1000       | False         | 527307     |
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

    ## define the groups
    # Each entry is a 4-tuple of lambda functions representing:
    # (name, order, gens, conditions)
    
    # groups which involve m,n
    chiral_1 = [
        # Table 4.1
        (
            lambda m,n: f"+-1/2[D{2*m}xD'{4*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_1, quat_j),
                         FDR(quat_j, quat_2D(1/(2*n)))
                         ),
            lambda m, n: m >= 2 and n >= 2,
        ),
        (
            lambda m,n: f"+-[D{2*m}xC{n}]",
            lambda m,n: 4*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_j, quat_1),
                         FDR(quat_1, quat_2D(1/n))
                         ),
            lambda m, n: m >= 2,
        ),
        (
            lambda m,n: f"+-1/2[D{2*m}xC{2*n}]",
            lambda m,n: 4*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_j, quat_2D(1/(2*n)))
                         ),
            lambda m, n: m >= 2 and n >= 2,
        ),
        (
            lambda m,n: f"+1/2[D{2*m}xC{2*n}]",
            lambda m,n: 2*m*n,
            lambda m,n: (FDR(quat_2D(1/m)*quat_minus1, quat_1),
                         FDR(quat_1, quat_2D(1/n)*quat_minus1),
                         FDR(quat_j, quat_2D(1/(2*n)))
                         ),
            lambda m,n: m % 2 == 1 and n % 2 == 1 and m >= 2 and n >= 2,
        ),
        (
            lambda m,n: f"+-1/2[D'{4*m}xC{2*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_j, quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_2D(1/(2*m)), quat_2D(1/(2*n)))
                         ),
            lambda m, n: m >= 2,
        ),

        # Mirrors of Table 4.1
        (
            lambda m,n: f"+-1/2[D'{4*m}xD{2*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_2D(1/m), quat_1),
                         FDR(quat_j, quat_1),
                         FDR(quat_2D(1/(2*m)), quat_j)
                         ),
            lambda m, n: n >= 2 and m >= 2,
        ),
        (
            lambda m,n: f"+-[C{m}xD{2*n}]",
            lambda m,n: 4*m*n,
            lambda m,n: (FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_1, quat_j),
                         FDR(quat_2D(1/m), quat_1)
                         ),
            lambda m, n: n >= 2,
        ),
        (
            lambda m,n: f"+-1/2[C{2*m}xD{2*n}]",
            lambda m,n: 4*m*n,
            lambda m,n: (FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_2D(1/m), quat_1),
                         FDR(quat_2D(1/(2*m)), quat_j)
                         ),
            lambda m, n: n >= 2 and m >= 2,
        ),
        (
            lambda m,n: f"+1/2[C{2*m}xD{2*n}]",
            lambda m,n: 2*m*n,
            lambda m,n: (FDR(quat_1, quat_2D(1/n)*quat_minus1),
                         FDR(quat_2D(1/m)*quat_minus1, quat_1),
                         FDR(quat_2D(1/(2*m)), quat_j)
                         ),
            lambda m,n: m % 2 == 1 and n % 2 == 1 and n >= 2 and m >= 2,
        ),
        (
            lambda m,n: f"+-1/2[C{2*m}xD'{4*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_1, quat_j),
                         FDR(quat_2D(1/m), quat_1),
                         FDR(quat_2D(1/(2*m)), quat_2D(1/(2*n)))
                         ),
            lambda m, n: n >= 2,
        ),

        # Table 4.2
        (
            lambda m,n: f"+-[D{2*m}xD{2*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_j, quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_1, quat_j)
                         ),
            lambda m,n: m >= 2 and n >= 2
        ),
        (
            lambda m,n: f"+-1/2[D'{4*m}xD'{4*n}]",
            lambda m,n: 16*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_j, quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_1, quat_j),
                         FDR(quat_2D(1/(2*m)), quat_2D(1/(2*n)))
                         ),
            lambda m,n:  m >= 2 and n >= 2,
        ),
        (
            lambda m,n: f"+-1/4[D{4*m}xD'{4*n}]",
            lambda m,n: 8*m*n,
            lambda m,n: (FDR(quat_2D(1/m), quat_1),
                         FDR(quat_1, quat_2D(1/n)),
                         FDR(quat_2D(1/(2*m)), quat_j),
                         FDR(quat_j, quat_2D(1/(2*n)))
                         ),
            lambda m,n:  m >= 2 and n >= 2,
        ),
        (
            lambda m,n: f"+1/4[D{4*m}xD'{4*n}]",
            lambda m,n: 4*m*n,
            lambda m,n: (FDR(quat_2D(1/m)*quat_minus1, quat_1),
                         FDR(quat_1, quat_2D(1/n)*quat_minus1),
                         FDR(quat_2D(1/(2*m)), quat_j),
                         FDR(quat_j, quat_2D(1/(2*n)))
                         ),
            lambda m,n: m % 2 == 1 and n % 2 == 1 and m >= 2 and n >= 2,
        ),
    ]

    ## groups which involve m,n,f,s
    chiral_2 = [
        ## Table 4.2
        (
            lambda m,n,f,s: f"+-1/{2*f}[D{2*m*f}xD{2*n*f}^{s}]",
            lambda m,n,f,s: 4*m*n*f,
            lambda m,n,f,s: (FDR(quat_2D(1/m), quat_1),
                             FDR(quat_1, quat_2D(1/n)),
                             FDR(quat_2D(1/(m*f)), quat_2D(s/(n*f))),
                             FDR(quat_j, quat_j)),
            lambda m,n,f,s: gcd(s,f) == 1 and 0 <= s <= f/2 and not(f == 1 and m == 1 and n == 1),
        ),
        (
            lambda m,n,f,s: f"+1/{2*f}[D{2*m*f}xD{2*n*f}^{s}]",
            lambda m,n,f,s: 2*m*n*f,
            lambda m,n,f,s: (FDR(quat_2D(1/m)*quat_minus1, quat_1),
                             FDR(quat_1, quat_2D(1/n)*quat_minus1),
                             FDR(quat_2D(1/(m*f)), quat_2D(s/(n*f))),
                             FDR(quat_j, quat_j)
                             ),
            lambda m,n,f,s: gcd(s,2*f) == 1 and m % 2 == 1 and n % 2 == 1 and 0 <= s <= f and not(f == 1 and m == 1 and n == 1),
        ),
        (
            lambda m,n,f,s: f"+-1/{f}[C{m*f}xC{n*f}^{s}]",
            lambda m,n,f,s: 2*m*n*f,
            lambda m,n,f,s: (FDR(quat_2D(1/m), quat_1),
                             FDR(quat_1, quat_2D(1/n)),
                             FDR(quat_2D(1/(m*f)), quat_2D(s/(n*f)))
                             ),
            lambda m,n,f,s: gcd(s,f) == 1 and 0 <= s <= f/2,
        ),
        (
            lambda m,n,f,s: f"+1/{f}[C{m*f}xC{n*f}^{s}]",
            lambda m,n,f,s: m*n*f,
            lambda m,n,f,s: (FDR(quat_2D(1/m)*quat_minus1, quat_1),
                             FDR(quat_1, quat_2D(1/n)*quat_minus1),
                             FDR(quat_2D(1/(m*f)), quat_2D(s/(n*f)))
                             ),
            lambda m,n,f,s: gcd(s,2*f) == 1 and m % 2== 1 and n % 2 == 1 and 0 <= s <= f,
        ),
    ]

    print("start generating chiral toroidal (CS) groups up to order", orderbound)

    ## chiral with 2 parameters m and n
    for typ in chiral_1:
        for m in IntegerRange(1, 1+orderbound):
            for n in IntegerRange(1, 1+orderbound):
                order = typ[1](m=m,n=n)
                if order > orderbound: break
                name = typ[0](m=m,n=n)
                gens = typ[2](m=m,n=n)
                conditions = typ[3]
                if not generate_full:
                    if conditions is not None and not conditions(m=m,n=n):
                        continue
                else:
                    if not name.startswith("+-"):
                        if is_even(m) or is_even(n):
                            continue
                gr = AxB_group(gens=gens, name=name)
                gr.latex = CS_name_to_latex(name)
                assert gr.order/2 == order
                insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
                num_groups += 1

    ## chiral with 4 parameters m, n, f, s
    for typ in chiral_2:
        for m in IntegerRange(1, 1+orderbound):
            for n in IntegerRange(1, 1+orderbound):
                for f in IntegerRange(1, 1+orderbound):
                    order = typ[1](m=m,n=n,f=f,s=None)
                    if order > orderbound: break
                    for s in IntegerRange(2*f): 
                        name = typ[0](m=m,n=n,f=f,s=s)
                        gens = typ[2](m=m,n=n,f=f,s=s)
                        conditions = typ[3]
                        if not generate_full:
                            if conditions is not None and not conditions(m=m,n=n,f=f,s=s): continue
                        else:
                            # Haploid groups require m and n being odd to be generated
                            if not name.startswith("+-"):
                                if is_even(m) or is_even(n) or gcd(s, 2*f) != 1:
                                    continue
                        gr = AxB_group(gens=gens, name=name)
                        gr.latex = CS_name_to_latex(name)
                        assert gr.order/2 == order
                        insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
                        num_groups += 1

    print("finished generating chiral toroidal groups (CS)")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")
    return

def generate_achiral_toroidal_CS(orderbound=130, generate_full=False,
                                 names_only=False, latex_only=False,
                                 verbose=True):
    """
    Generate achiral toroidal groups (CS) and insert them in `catalog`.

    INPUT:
    - `orderbound` -- int: Order bound on the generated group.
    - `generate_full` -- boolean: Whether to include duplications or not.
    - `names_only` -- boolean: If True, inserted values in `catalogs` are
       names of the groups.
    - `latex_only` -- boolean: If True, inserted values in `catalogs` are 
       latex names of the groups. This option overrides `names_only`.
    - `verbose` -- boolean: More information.

    WARNING:
    The option `generate_full` produces duplications!

    REMARK:
    `catalog` should be defined before running this function. Otherwise,
    you get a `NameError`.

    RUNNING TIME:
    orderbound | generate_full | num_groups | time
    ---------- | ------------- | ---------- | ------------------
    130        | True          |            |                    
    130        | False         |            |                   
    500        | True          |            |                    
    500        | False         |            |                    
    1000       | False         |            |
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

    # groups which involve n
    achiral_1 = [
        # Table 4.3
        (
            lambda n: f"+-[D{2*n}xD{2*n}].2",
            lambda n: 16*n^2,
            lambda n: (FDR(quat_2D(1/n), quat_1),
                       FDR(quat_j, quat_1),
                       FDR(quat_1, quat_2D(1/n)),
                       FDR(quat_1, quat_j),
                       FDR(quat_1, quat_1, True)
                       ),
            None,
            None,
        ),
        (
            lambda n: f"+-1/2[D'{4*n}xD'{4*n}].2",
            lambda n: 32*n^2,
            lambda n: (FDR(quat_2D(1/n), quat_1),
                       FDR(quat_j, quat_1),
                       FDR(quat_1, quat_2D(1/n)),
                       FDR(quat_1, quat_j),
                       FDR(quat_2D(1/(2*n)), quat_2D(1/(2*n))),
                       FDR(quat_1, quat_1, True)
                       ),
            None,
            None,
        ),
        (
            lambda n: f"+-1/2[D'{4*n}xD'{4*n}].-2",
            lambda n: 32*n^2,
            lambda n: (FDR(quat_2D(1/n), quat_1),
                       FDR(quat_j, quat_1),
                       FDR(quat_1, quat_2D(1/n)),
                       FDR(quat_1, quat_j),
                       FDR(quat_2D(1/(2*n)), quat_2D(1/(2*n))),
                       FDR(quat_1, quat_2D(1/(2*n)), True)
                       ),
            None,
            None,
        ),
        (
            lambda n: f"+-1/4[D{4*n}xD'{4*n}].2",
            lambda n: 16*n^2,
            lambda n: (FDR(quat_2D(1/n), quat_1),
                       FDR(quat_1, quat_2D(1/n)),
                       FDR(quat_2D(1/(2*n)), quat_j),
                       FDR(quat_j, quat_2D(1/(2*n))),
                       FDR(quat_1, quat_1, True)
                       ),
            None,
            None,
        ),
        (
            lambda n: f"+1/4[D{4*n}xD'{4*n}].2_3",
            lambda n: 8*n^2,
            lambda n: (FDR(quat_2D(1/n)*quat_minus1, quat_1),
                       FDR(quat_1, quat_2D(1/n)*quat_minus1),
                       FDR(quat_2D(1/(2*n)), quat_j),
                       FDR(quat_j, quat_2D(1/(2*n))),
                       FDR(quat_1, quat_1, True)
                       ),
            lambda n: n % 2 == 1,
            None,
        ),
        (
            lambda n: f"+1/4[D{4*n}xD'{4*n}].2_1",
            lambda n: 8*n^2,
            lambda n: (FDR(quat_2D(1/n)*quat_minus1, quat_1),
                       FDR(quat_1, quat_2D(1/n)*quat_minus1),
                       FDR(quat_2D(1/(2*n)), quat_j),
                       FDR(quat_j, quat_2D(1/(2*n))),
                       FDR(quat_minus1, quat_1, True)),
            lambda n: n % 2 == 1,
            None,
        ),
    ]

    ## achiral with 1 parameter
    for typ in achiral_1:
        for n in IntegerRange(1, 1+orderbound):
            order = typ[1](n=n)
            if order > orderbound: break

            gen_conditions = typ[3]
            dup_conditions = typ[4]
            if gen_conditions is not None and not gen_conditions(n=n): 
                continue
            if not generate_full:
                if dup_conditions is not None and not dup_conditions(n=n):
                    continue

            gens = typ[2](n=n)
            name = typ[0](n=n)
            gr = AxB_group(gens=gens, name=name)

            gr.latex = CS_name_to_latex(name)
            assert gr.order/2 == order
            insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
            num_groups += 1

    # groups which involve n,f,s,alpha,beta
    achiral_2 = [
        # Table 4.3
        (
            lambda n,f,s,alpha,beta: f"+-1/{2*f}[D{2*n*f}xD{2*n*f}^{s}].2^({alpha},{beta})",
            lambda n,f,s,alpha,beta: 8*n^2*f,
            lambda n,f,s,alpha,beta: (FDR(quat_2D(1/n), quat_1),
                                      FDR(quat_1, quat_2D(1/n)),
                                      FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                                      FDR(quat_j, quat_j),
                                      FDR(quat_2D(alpha/(2*n*f)), quat_2D((alpha*s+beta*f)/(2*n*f)), True)
                                      ),
            lambda n,f,s,alpha,beta: gcd(s, f) == 1 and\
                                     ((s^2-1)/f).is_integer() and\
                                     (alpha*((s^2-1)/f) + beta*f) % 2 == 0 and\
                                     (alpha, beta) in [(0,0), (0,1), (1,0), (1,1)],
            lambda n,f,s,alpha,beta: 0 <= s <= f/2,
        ),
        (
            lambda n,f,s,alpha,beta: f"+1/{2*f}[D{2*n*f}xD{2*n*f}^{s}].2^({alpha},{beta})",
            lambda n,f,s,alpha,beta: 4*n^2*f,
            lambda n,f,s,alpha,beta: (FDR(quat_2D(1/n)*quat_minus1, quat_1),
                                      FDR(quat_1, quat_2D(1/n)*quat_minus1),
                                      FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                                      FDR(quat_j, quat_j),
                                      FDR(quat_2D(alpha/(2*n*f)), quat_2D((alpha*s+beta*f)/(2*n*f)), True)
                                      ),
            lambda n,f,s,alpha,beta: gcd(s,2*f) == 1 and is_odd(n) and\
                                     ((s^2-1)/f).is_integer() and\
                                     (alpha*((s^2-1)/f)) % 4 == 0 and\
                                     is_even((s^2-1)/f) and\
                                     beta % 2 == 0 and\
                                     (alpha, beta) in [(0,0), (0,2), (1,0), (1,2)],
            lambda n,f,s,alpha,beta: 0 <= s <= f,
        ),
    ]

    ## achiral with 5 parameters
    for typ in achiral_2:
        for n in IntegerRange(1, 1+orderbound):
            for f in IntegerRange(1, 1+orderbound):
                order = typ[1](n=n, f=f, s=None, alpha=None, beta=None)
                if order > orderbound: break
                for s in IntegerRange(2*f):
                    for alpha in IntegerRange(2):
                        for beta in IntegerRange(3):
                            gen_conditions = typ[3]
                            dup_conditions = typ[4]
                            if gen_conditions is not None and not gen_conditions(n=n, f=f, s=s, alpha=alpha, beta=beta): 
                                continue
                            if not generate_full:
                                if dup_conditions is not None and not dup_conditions(n=n, f=f, s=s, alpha=alpha, beta=beta): continue
                            gens = typ[2](n=n, f=f, s=s, alpha=alpha, beta=beta)
                            name = typ[0](n=n, f=f, s=s, alpha=alpha, beta=beta)
                            gr = AxB_group(gens=gens, name=name)
                            gr.latex = CS_name_to_latex(name)
                            assert gr.order/2 == order, (n, f, s, alpha, beta)
                            insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
                            num_groups += 1

    # groups which involve n,f,s
    achiral_3 = [
        # Table 4.3
        (
            lambda n,f,s: f"+-1/{2*f}[D{2*n*f}xD{2*n*f}^{s}].-2",
            lambda n,f,s: 8*n^2*f,
            lambda n,f,s: (FDR(quat_2D(1/n), quat_1),
                           FDR(quat_1, quat_2D(1/n)),
                           FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                           FDR(quat_j, quat_j),
                           FDR(quat_1, quat_j, True)
                           ),
            lambda n,f,s: ((s^2+1)/f).is_integer() and\
                          gcd(f,s) == 1,
            lambda n,f,s: 0 <= s <= f/2,
        ),
        (
            lambda n,f,s: f"+1/{2*f}[D{2*n*f}xD{2*n*f}^{s}].-2",
            lambda n,f,s: 4*n^2*f,
            lambda n,f,s: (FDR(quat_2D(1/n)*quat_minus1, quat_1), 
                           FDR(quat_1, quat_2D(1/n)*quat_minus1),
                           FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                           FDR(quat_j, quat_j),
                           FDR(quat_1, quat_j, True)
                           ),
            lambda n,f,s: ((s^2+1)/f).is_integer() and\
                          is_even((s^2+1)/f) and\
                          n%2 == 1 and\
                          gcd(s,2*f) == 1,
            lambda n,f,s: 0 <= s <= f,
        ),
    ]


    ## achiral with 3 parameters
    for typ in achiral_3:
        for n in IntegerRange(1, 1+orderbound):
            for f in IntegerRange(1, 1+orderbound):
                order = typ[1](n=n, f=f, s=None)
                if order > orderbound: break
                for s in IntegerRange(2*f):
                    gen_conditions = typ[3]
                    dup_conditions = typ[4]
                    if gen_conditions is not None and not gen_conditions(n=n, f=f, s=s): 
                        continue
                    if not generate_full:
                        if dup_conditions is not None and not dup_conditions(n=n, f=f, s=s):
                            continue
                    gens = typ[2](n=n, f=f, s=s)
                    name = typ[0](n=n, f=f, s=s)
                    gr = AxB_group(gens=gens, name=name)
                    gr.latex = CS_name_to_latex(name)
                    assert gr.order/2 == order
                    insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
                    num_groups += 1


    # groups which involve n,f,s,gamma
    achiral_4 = [   
        # Table 4.3
            (
                lambda n,f,s,gamma: f"+-1/{f}[C{n*f}xC{n*f}^{s}].2^({gamma})",
                lambda n,f,s,gamma: 4*n^2*f,
                lambda n,f,s,gamma: (FDR(quat_2D(1/n), quat_1),
                                     FDR(quat_1, quat_2D(1/n)),
                                     FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                                     FDR(quat_1, quat_2D(gamma*gcd(f,s+1)/(2*n*f)), True)
                                     ),
                # modified condition
                lambda n,f,s,gamma: gcd(s,f) == 1 and\
                                    ((s^2-1)/f).is_integer() and\
                                    (gcd((s^2-1)/f, s-1)*gamma) % 2 == 0 and\
                                    (gcd(f, s+1)*gamma) % 2 == 0,
                lambda n,f,s,gamma: 0 <= s < f,
            ),
            (
                lambda n,f,s,gamma: f"+1/{f}[C{n*f}xC{n*f}^{s}].2^({gamma})",
                lambda n,f,s,gamma: 2*n^2*f,
                lambda n,f,s,gamma: (FDR(quat_2D(1/n)*quat_minus1, quat_1),
                                     FDR(quat_1, quat_2D(1/n)*quat_minus1),
                                     FDR(quat_2D(1/(n*f)), quat_2D(s/(n*f))),
                                     FDR(quat_1, quat_2D(gamma*gcd(f,s+1)/(2*n*f)), True)
                                     ),
                # modified condition
                lambda n,f,s,gamma: gcd(s,2*f) == 1 and\
                                    n % 2 == 1 and\
                                    ((s^2-1)/f).is_integer() and\
                                    (gcd((s^2-1)/f, s-1)*gamma) % 4 == 0 and\
                                    (gcd(f, s+1)*gamma) % 2 == 0 and\
                                    is_even((s^2-1)/f),
                lambda n,f,s,gamma: 0 <= s < 2*f,
            ),
            
    ]

    ## achiral with 4 parameters
    first = True        # these two groups have different gamma limit
    for typ in achiral_4:
        for n in IntegerRange(1, 1+orderbound):
            for f in IntegerRange(1, 1+orderbound):
                order = typ[1](n=n, f=f, s=None, gamma=None)
                if order > orderbound: break
                for s in IntegerRange(2*f):
                    if first:
                        gamma_values = [0,1]
                    else:
                        d = ZZ(gcd(2*f,s+1)/gcd(f,s+1))
                        gamma_limit = 2*d
                        gamma_values = [0,d]
                    for gamma in gamma_values:

                        gen_conditions = typ[3]
                        dup_conditions = typ[4]
                        if gen_conditions is not None and not gen_conditions(n=n, f=f, s=s, gamma=gamma): 
                            continue
                        if not generate_full:
                            if dup_conditions is not None and not dup_conditions(n=n, f=f, s=s, gamma=gamma):
                                continue
                        gens = typ[2](n=n, f=f, s=s, gamma=gamma)
                        name = typ[0](n=n, f=f, s=s, gamma=gamma)
                        gr = AxB_group(gens=gens, name=name)
                        gr.latex = CS_name_to_latex(name)
                        assert gr.order/2 == order
                        insert_into_catalog(gr, names_only=names_only, latex_only=latex_only, verbose=verbose)
                        num_groups += 1
        first = False

    print("finished generating achiral toroidal groups (CS)")
    print(f"number of generated groups: {num_groups}")
    print(f"number of unique groups in catalog: {len(catalog)}")
    print(f"potential duplications {sum(map(len,catalog.values()))-len(catalog)}")
    print(f"elapsed time: {time()-start_time}s")

    return


def make_toroidal_cross_reference(orderbound, ordered=False):
    """
    Create a cross reference between CS names and our names for tubular
    groups. Store output in two files:
    * `logs/toroidal_chiral_cross_reference.log`
    * `logs/toroidal_achiral_cross_reference.log`
    """
    global catalog
    
    func_options =  {"orderbound": orderbound, 
                     "generate_full": False,
                     "names_only": False,
                     "latex_only": False,
                     "verbose": False}

    for func, func_CS, name in [(generate_chiral_toroidal, generate_chiral_toroidal_CS, "chiral"),
                                (generate_achiral_toroidal, generate_achiral_toroidal_CS, "achiral")]:
        catalog = {}
        func_CS(**func_options)
        catalog_CS = catalog.copy()

        catalog = {}
        func(**func_options)

        for cat in catalog:
            if cat not in catalog_CS:
                raise f"Some of our {name} groups is missing from CS"
            catalog[cat] += catalog_CS[cat]
        assert len(catalog) == len(catalog_CS), f"Some of CS {name} groups is missing from ours"


        outfile = open(f"logs/toroidal_{name}_cross_reference.log", "w")
        outfile.write(f"# groups up to order: {orderbound}\n")
        outfile.write(f"ord : our group       == CS group(s)\n")

        if ordered: 
            catalog_values = sorted(catalog.values(), key= lambda gs: gs[0].order)
        else:
            catalog_values = catalog.values()
        for gs in catalog_values:
            _temp = outfile.write(f"{int(gs[0].order/2):<4}: {gs[0].name:<14} ")
            for g in gs[1:]:
                _temp = outfile.write(f" == {g.name}")
            _temp = outfile.write("\n")
        outfile.close()

    return
