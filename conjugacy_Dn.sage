#!/usr/bin/env python # to trigger proper syntax highlighting

# TO BE CLEANED
# CAREFULL: it overwrites e (to be fixed)

def classify_Dn_group(g):
    """return "C",n,None or "D",2n,offset
    g is an iterator over the l or r components of the
    orientation-preserving elements.

    For Dn, the offset is the smallest alpha such that quat_2D(alpha,1) is in g
    """
    g = set(g)
    typ = "C"
    offset = 99
    for s in g:
        if s.p:
            typ = "D"
            offset = min(offset,s.alpha)
    if typ=="C":
        assert len(g)%2==0
        return "C",len(g)/2,None
    assert len(g)%4==0
    return "D",len(g)/2,offset

def group_table(g):
    g_chiral = [e for e in g if not e.star]
    achiral = len(g_chiral)<len(g)
    d1,n1,_= classify_Dn_group(e.l for e in g_chiral)
    d2,n2,_= classify_Dn_group(e.r for e in g_chiral)
    mod1 = n1 if d1=="D" else 2*n1 # elements are multiples of 2*pi/mod1
    mod2 = n2 if d2=="D" else 2*n2

    chiraltext = ""
    if not achiral:
        chiraltext = " (chiral)"
    
    present = {(i,j,star):' ' for i in range(2*n1) for j in range(2*n2) for star in (False,True)}
    for e in g:
        try:
            present[el_clas(e.l,mod1),el_clas(e.r,mod2),e.star] = "X"
        except AssertionError:
            if e.star:
                if achiral:
                    print("group element not good",e,n1,n2)
                chiraltext = " (only the chiral part)"
                achiral = False
            else:
                print("group element not good",e,n1,n2)
                print(" ".join(sorted(str(x) for x in g)))
                return
    top= "".join(str(i%10) for i in range(mod2))
    mid= mod2*"-"
    if d2=="D":
        top += "|"+top
        mid += "+"+mid
    if achiral:
        top += " *"+top
        mid += "  "+mid
 
    result = ["  "+top+chiraltext+f" {d1}{n1}x{d2}{n2}"]
    for i in range(2*mod1):
        if i==mod1:
            if d1=="C":
                break
            result.append("  "+mid)
        line = ""
        for star in False,True:
            line += " "
            line += (str((i%(mod1))%10)
                     +"".join(present[i,j,star] for j in range(mod2))+
                     ("|"+"".join(present[i,j,star] for j in range(n2,2*n2)) if d2=="D" else ""))      
            if not achiral:
                break
        result.append(line)
    for line in result:
        print (line)

def el_clas(e,n):
    "index in the table"
    res =e.alpha*n/2
    assert res.is_integer()
    assert 0<=res<n
    if e.p:
        return n+res
    return res


WHATHAPPENED = set()
STATISTICS = Counter()
STATISTICS["OLD METHOD "]=0

def perform_conjugates(what, conj, g, conju, name, verbose):
    if verbose:
        print(what, name, "conjugation with", conj)
    g = [AxB_group.conjugate_element(e,*conj) for e in g]
    if verbose>2:
        print(" ".join(sorted(map(str,g))))
    if verbose>1:
        group_table(g)
    conju.append(conj)
    return g,conju

def normalize_orientation_reversing_part(g,conju,name="",verbose=2):
    """ensure that all elements l,r of [l,r] or *[l,r] belong to D_2n or C_n
    (Step 1 of the write-up on extending element for achiral groups:
    conjugation with [1,a], for some element *[a,b] in g

    prefer an element a with a.p==False"""
    achiral = False
    cyclic_found = False
    for e in g:
        if e.star:
            if e.l.p:
                if not cyclic_found:
                    a = e.l
            else:
                a = e.l
                cyclic_found = True
            achiral = True
            if e.l==quat_1:
                return g,conju,achiral # no conjugation necessary
    if achiral:
        g,conju = perform_conjugates("normalize orientation-reversing part",
                                     (quat_1,a),
                                     g,conju,name,verbose)
    return g,conju,achiral
    
def normalize_dihedral_part(g,conju,typ,name="",verbose=2):

    if typ[0]=="D" and typ[2]!=0:
        start1 = quat_2D(typ[2]/2,0)
        return perform_conjugates("rotating",
                              (start1,start1),
                              g,conju,name,verbose)
    return g,conju

def normalize_dihedral_part_chiral(g,conju,name="",verbose=2):
    offset_l = offset_r = 99
    for e in g:
        if e.l.p:
            offset_l = min(offset_l,e.l.alpha)
        if e.r.p:
            offset_r = min(offset_r,e.r.alpha)
    if offset_l == offset_r == 99:
        return g,conju
    if offset_l == 99: offset_l = 0 
    if offset_r == 99: offset_r = 0 

    if offset_l == offset_r == 0:
        return g,conju
    else:
        return perform_conjugates("rotating",
                              (quat_2D(offset_l/2),quat_2D(offset_r/2)),
                              g,conju,name,verbose)

def change_C2_to_D2(g,conju,name="",verbose=2):
    typl = classify_Dn_group(e.l for e in g if not e.star)
    typr = classify_Dn_group(e.r for e in g if not e.star)
    if typl[:2]==("C",2) or typr[:2]==("C",2):
        conl = ("S",0) if typl[:2]==("C",2) else ("S",-1)
        conr = ("S",0) if typr[:2]==("C",2) else ("S",-1)
        g,conju = perform_conjugates("converted C2->D2 group:",
                                  (conl,conr),
                                  g,conju,name,verbose)
        g,conju,achiral = normalize_orientation_reversing_part(g,conju, name=name,verbose=verbose)
          # is this really necessary?
    return g,conju

def fix_cyclic_part(g_reference, g,conju,name="",verbose=2):
    cyclic1 = set(e for e in g_reference if not e.star and not e.l.p and not e.r.p)
    cyclic2 = set(e for e in g           if not e.star and not e.l.p and not e.r.p)
    if cyclic1 != cyclic2:
        what = "orientation reversal cyclic part G2:"
        STATISTICS[what]+=1
        return perform_conjugates(what,
                                  (quat_j,quat_1),
                                  g,conju,name,verbose)
        # this might be redundant in case of a D4-subgroup.
    return g,conju

class Identity: # replacement for tr[..] when not a D4-group
    def __getitem__(self,k):
        return(k)
id = Identity()   

good_pairs = Counter()
good_pairs2 = Counter()

def try_D4_special_chiral_first(g1,g2,isL,isR):
    ## special treatment in advance:
    tr,con_alt = transform_D4[1]
    g1_set = set(g1)
    if isR:
        g2b = [FDR(e.l,tr[e.r]) for e in g2] # the most frequent version
        if all(e in g1_set for e in g2b):
            return (quat_1,con_alt),g2b
        if isL:
            g2b = [FDR(tr[e.l],tr[e.r]) for e in g2]
            if all(e in g1_set for e in g2b):
                return (con_alt,con_alt),g2b
    if isL:
        g2b = [FDR(tr[e.l],e.r) for e in g2]
        if all(e in g1_set for e in g2b):
            return (con_alt,quat_1),g2b
    return None,g2

STATISTICS["special_chiral_not_enough"]=0

def try_D4_special_chiral(g1,g2,conju1,conju2,typl,typr, verbose):
    #OK, g1q, g2q,con1,con2 = try_D4_chiral_special(...

    ## special treatment in advance:
    conj2,g2 = try_D4_special_chiral_first(g1,g2, is_subgroup_of_D4(typl),is_subgroup_of_D4(typr))
    if conj2 is not None:
        good_pairs[str(conj2)] += 1
        conju2.append(conj2)
        return True,g1,g2,conju1,conju2

    STATISTICS["special_chiral_not_enough"]+=1
            
    ## standard treatment. (This should never be used)
    if is_subgroup_of_D4(typl):
        loop_left = transform_D4
    else:
        loop_left = (id,quat_1), # the identity
    if is_subgroup_of_D4(typr):
        loop_right = transform_D4
    else:
        loop_right = (id,quat_1),
    g1_set = set((e.l,e.r) for e in g1)
    for tr2l,con2l in loop_left:
        for tr2r,con2r in loop_right:
            # the redundant pair (id,id) is always checked first
            for e in g2:
                if (tr2l[e.l],tr2r[e.r]) not in g1_set:
                    break
            else: # found
                conj2 = con2l,con2r
                if verbose:
                    print ("special D4: G2 conjugated by", conj2)
                g2 = [FDR(tr2l[e.l],tr2r[e.r]) for e in g2]
                assert set(g1) == set(g2),"sets should be equal"
                if verbose>1:
                    print("G1:")
                    group_table(g1)                    
                    print("G2:")
                    group_table(g2)                    
                good_pairs[str(conj2)] += 1
                conju2.append(conj2)
                return True,g1,g2,conju1,conju2
    return False,g1,g2,conju1,conju2

def try_D4_special_achiral(g1,g2,conju1,conju2,name1,name2, verbose):
    g1_set = set(g1)
    for it,(tr2,conj2) in enumerate(transform_D4xD4):
        if all(tr2[e] in g1_set for e in g2):
            g2 = [tr2[e] for e in g2]
            if verbose:
                print("special D4 achiral: G2 conjugated by",conj2)
            conju2.append(conj2)
            good_pairs2[it] += 1
            if it>=8:
                WHATHAPPENED.add(f"special_achiral_not_enough {it} {name1} {name2}")
            ## This was only for safety checking:
            # g2 = [AxB_group.conjugate_element(e.to_exact_quaternions(), *conj2) for e in g2]
            # g1 = [e.to_exact_quaternions() for e in g1] # to be on equal footing
            return True,g1,g2,conju1,conju2
    return False,g1,g2,conju1,conju2

DIAG = Counter()
def normalize_noncyclic_part_achiral(g, conju, typ, name="",verbose=2):
    "within the achiral part"
    alpha = alpha_r = alpha_l = alpha_lr = False # infinity
    alpha_r0 = 99 # infinity
    for e in g:
        if not e.star: continue
        if not e.r.p and not e.l.p :
            alpha = True
        if not e.r.p and e.l.p :
            alpha_l = True
        if e.r.p and not e.l.p :
            alpha_r = True
        if e.l.p and e.r.p :
            alpha_lr = True
        if e.l==quat_1 and e.r.p :
            alpha_r0 = min(alpha_r0, e.r.alpha)
    #print(f"{alpha=} {alpha_l=} {alpha_r=} {alpha_lr=} {alpha_r0=}",name)
    if alpha_l == alpha_r ==alpha == alpha_lr == False:
        return (g,conju) # chiral group
    if not alpha and not alpha_lr and alpha_r0 < 99:
        # group L: easy
        con_l = con_r = alpha_r0/2
    else:
        d,n,_ = typ
        mod1 = n if d=="D" else 2*n # elements are multiples of 2/mod1 (actually, 2*pi/mod1)
        if alpha:
            # group | or +
            diag = 2*mod1 # infinity
            for e in g:
                if not e.star: continue
                if not e.l.p and not e.r.p:
                    diag_new = int(mod1*(e.l*e.r).alpha) # should divide by 2 to get the old diag
                    if diag_new < diag:
                        diag = diag_new
                        alpha0 = e.l.alpha
                        #print(f"{alpha0=} {diag=} {e=}")
                    elif diag_new == diag:
                        alpha0 = min(alpha0,e.l.alpha)
            DIAG["0",diag//2] += 1            
        else:
            alpha0 = 0

        if alpha_lr:
            # group - or +
            alpha1 = 99           
            for diag in range(mod1):
                # look for first occupied diagonal l*r^-1=diag
                for e in g:
                    if not e.star: continue
                    if e.l.p and e.r.p and e.l*e.r.conjugate() == quat_2D(2*diag/mod1):
                        alpha1 = min(e.l.alpha, alpha1)
                        #print(f"{alpha1=} {diag=} {e=}")
                if alpha1 != 99: 
                    DIAG["1",diag] += 1
                    break
        else:
            alpha1 = 0
        #print(f"**OR(2) {alpha0 = } {alpha1 = }",name)
            
        # if alpha0 is not None and alpha1 is not None:
        #   group "+":
        #   solve
        #       con_l - con_r == - alpha0
        #      -con_l - con_r == - alpha1
        #   in order to reach FDR(quat_2D(0,*),*,True)            
        # elif alpha0 is not None:
        #     group "|", type "C":
        # elif alpha1 is not None:
        #     group "-", type "C":
            
        con_l = (-alpha0+alpha1)/2
        con_r =  (alpha0+alpha1)/2

    what = "normalize non-cyclic part (achiral, 2)"
        
    if con_l == con_r == 0:
        if verbose: print(what,name,"nothing to do")
        return g,conju
    return perform_conjugates(what,
                              (quat_2D(con_l),quat_2D(con_r)),
                              g,conju,name,verbose)

def normalize_noncyclic_part_chiral(g, conju, name="",verbose=2):
    "within the chiral part"
    alpha_r = alpha_l = alpha_lr = 99 # infinity
    beta_r = beta_l = beta_lr = 99 # infinity
    for e in g: # find the lowest occupied line in the table
        if e.star: continue
        if not e.r.p and e.l.p :
            beta_l = min(beta_l, e.r.alpha)
        elif not e.l.p and e.r.p :
            beta_r = min(beta_r, e.l.alpha)
        elif e.l.p and e.r.p :
            beta_lr = min(beta_lr, e.l.alpha)
            
    for e in g: # find in the lowest occupied line the leftmost occupied entry
        if e.star: continue
        if not e.r.p and e.r.alpha==beta_l and e.l.p :
            alpha_l = min(alpha_l, e.l.alpha)
        if not e.l.p and e.l.alpha==beta_r and e.r.p :
            alpha_r = min(alpha_r, e.r.alpha)
        if e.l.p and e.l.alpha==beta_lr and e.r.p :
            alpha_lr = min(alpha_lr, e.r.alpha)

    #print(f"{alpha_r = } {alpha_l = } {alpha_lr = }")
    #print(f"{beta_r = } {beta_l = } {beta_lr = }")

    con_l = con_r = 0
    if alpha_l != 99 or alpha_r != 99:
        if alpha_l != 99:
            con_l = alpha_l/2
        if alpha_r != 99:
            con_r = alpha_r/2
    elif alpha_lr != 99:
        con_l,con_r = (0,alpha_lr/2)
        #con_l,con_r = (alpha_lr/4,alpha_lr/4) # ?? better for achiral groups

    what = "normalize non-cyclic part (chiral)"
    if con_l==con_r==0:
        if verbose>2: print(what,name, "nothing to do")
        return g,conju
    
    g,conju = perform_conjugates(what,
                              (quat_2D(con_l),quat_2D(con_r)),
                              g,conju,name,verbose)
    g,conju,achiral = normalize_orientation_reversing_part(g,conju, name=name,verbose=verbose)
        # seems to work somewhat by chance
    return g,conju

def try_diagonal_move(g, conju, name="",verbose=2): # for G*^gS_2
    # bring *[l,r] with l.p==r.p==False into the corner
    for e in g:
        if e.star and not e.l.p and not e.r.p:
            if e.l.alpha==e.r.alpha:
                con_l,con_r = e.l.alpha, 0
            break
    else:
        return g,conju
    what = 'try "diagonal move" for G*gS'
    WHATHAPPENED.add("diagonal-move: "+g.name)
    STATISTICS[what]+=1
    return perform_conjugates(what,
                              (quat_2D(con_l),quat_2D(con_r)),
                              g,conju,name,verbose)
    
        
def is_subgroup_of_D4(typ):
    d,n,_= typ
    return (d,n) in (("D",4),("D",2),("C",2),("C",1))

def find_conjugation_2D(g1,g2,name1="gr1",name2="gr2",verbose=0):

    assert len(g1)==len(g2)

#    achiral  = any(e.star for e in g1)
#    achiral2 = any(e.star for e in g2)
#    assert achiral==achiral2

    con1 = []
    con2 = []
    
    if 0:#achiral:
        g1,con1 = normalize_noncyclic_part(g1,con1,name=name1,verbose=verbose)
        g2,con2 = normalize_noncyclic_part(g2,con2,name=name2,verbose=verbose)

    
    g1,con1,achiral  = normalize_orientation_reversing_part(g1,con1,name=name1,verbose=verbose)
    g2,con2,achiral2 = normalize_orientation_reversing_part(g2,con2,name=name2,verbose=verbose)
    assert achiral==achiral2

    typ1 = classify_Dn_group(e.l for e in g1 if not e.star)
    typ2 = classify_Dn_group(e.l for e in g2 if not e.star)
    if verbose:
        print("Types",typ1,typ2)
        if verbose>2:
            print(name1," ".join(sorted(map(str,g1))))
            print(name2," ".join(sorted(map(str,g2))))
    if achiral:
        g1,con1 = normalize_dihedral_part(g1,con1,typ1,name=name1,verbose=verbose)
        g2,con2 = normalize_dihedral_part(g2,con2,typ2,name=name2,verbose=verbose)
    else:
        g1,con1 = normalize_dihedral_part_chiral(g1,con1,name=name1,verbose=verbose)
        g2,con2 = normalize_dihedral_part_chiral(g2,con2,name=name2,verbose=verbose)
    if verbose>1:
        print("group 1",name1, "after normalization")
        group_table(g1)
        if verbose>2:
            print(" ".join(sorted(map(str,g1))))
            print(",".join(sorted(map(repr,g1))))
        print("group 2",name2, "after normalization")
        group_table(g2)  
        if verbose>2:
            print(" ".join(sorted(map(str,g2))))
            print(",".join(sorted(map(repr,g2))))
        
    g1,con1 = change_C2_to_D2(g1,con1,name=name1,verbose=verbose)
    g2,con2 = change_C2_to_D2(g2,con2,name=name2,verbose=verbose)
    
    typ1 = classify_Dn_group(e.l for e in g1 if not e.star)
    typ2 = classify_Dn_group(e.l for e in g2 if not e.star)
    if verbose:
        print("Types",typ1,typ2)
    if typ1[0]!=typ2[0]:
            print("Types",typ1,typ2)
            print(typ1,typ2,"different groups. not implemented", name1, name2)
            raise NotImplementedError
            return
    g2,con2 = fix_cyclic_part(g1, g2,con2,name=name2,verbose=verbose)

    if set(g1)==set(g2): return True,con1,con2

    cyclic1 = set(e for e in g1 if not e.star and not e.l.p and not e.r.p)
    cyclic2 = set(e for e in g2 if not e.star and not e.l.p and not e.r.p)
    if cyclic1 != cyclic2:
        if achiral and (
                is_subgroup_of_D4(typ1) and is_subgroup_of_D4(typ2)): # subgroup of D4xD4
            if verbose>2:
                print(" ".join(sorted(map(str,g1))))
                print(",".join(sorted(map(repr,g1))))
                print(" ".join(sorted(map(str,g2))))
                print(",".join(sorted(map(repr,g2))))
            # try the prerecorded conjugacies
            OK, g1q, g2q,con1,con2 = try_D4_special_achiral(g1,g2, con1,con2,name1,name2,verbose=verbose)
            if OK:
                assert set(g1q) == set(g2q)
                return True,con1,con2
        if achiral:
            print("not found give up")
            return False, "not found give up"

    if 1: #not achiral:
        g1,con1 = normalize_noncyclic_part_chiral(g1,con1,name=name1,verbose=verbose)
        g2,con2 = normalize_noncyclic_part_chiral(g2,con2,name=name2,verbose=verbose)

    if set(g1)==set(g2): return True,con1,con2
    
    if achiral:
        g1,con1 = normalize_noncyclic_part_achiral(g1, con1, typ1, name=name1, verbose=verbose)
        g2,con2 = normalize_noncyclic_part_achiral(g2, con2, typ2, name=name2, verbose=verbose)

    if set(g1)==set(g2): return True, con1, con2

    if achiral:
        if is_subgroup_of_D4(typ1) and is_subgroup_of_D4(typ2): # subgroup of D4xD4
            # try the prerecorded conjugacies
            OK, g1q, g2q,con1,con2 = try_D4_special_achiral(g1,g2, con1,con2,name1,name2,verbose=verbose)
            if OK:
                assert set(g1q) == set(g2q)
                return True,con1,con2

        con20 = con2.copy()
        g2b,con2 = perform_conjugates("try (1,j)",
                                     (quat_1,quat_j),
                                     g2,con2,name2,verbose)
        if set(g1)==set(g2b):
            what = 'conjugate (1,j) achiral'
            WHATHAPPENED.add(what + ": "+name1)
            STATISTICS[what]+=1


            return True, con1, con2
        con2 = con20
        
        g1,con1 = try_diagonal_move(g1,con1,name=name1,verbose=verbose)
        g2,con2 = try_diagonal_move(g2,con2,name=name2,verbose=verbose)
        if set(g1)==set(g2b): return True, con1, con2
        
        print("achiral, give up")
        return False,
    else:        
        typ1l = classify_Dn_group(e.l for e in g1 if not e.star)
        typ1r = classify_Dn_group(e.r for e in g1 if not e.star)
        typ2l = classify_Dn_group(e.l for e in g2 if not e.star)
        typ2r = classify_Dn_group(e.r for e in g2 if not e.star)
        if verbose:
            print (f"G1 = {typ1l}x{typ1r}, G2 = {typ2l}x{typ2r}")
        if is_subgroup_of_D4(typ1l) or is_subgroup_of_D4(typ1r):
            if verbose>2:
                print("G1")
                print(" ".join(sorted(map(str,g1))))
                print("G2")
                print(" ".join(sorted(map(str,g2))))
            assert typ1l==typ2l and typ1r==typ2r
            # try the special conjugacies
            OK, g1q, g2q,con1,con2 = try_D4_special_chiral(g1, g2,con1,con2,typ1l,typ1r,verbose=verbose)
            if OK:
                return True,con1,con2
        
        print ('in "2 .',typ1,"x",typ2,'"', name1,name2,'   1xorder =',len(g1)/2)
        if 1:   
                print("G1")
                print(" ".join(sorted(map(str,g1))))
                print(name1, "=[")
                print(",".join(sorted(map(repr,g1))),"]")
                print("G2")
                print(" ".join(sorted(map(str,g2))))
                print(name2, "=[")
                print(",".join(sorted(map(repr,g2))),"]")
        return False,"not found give up"

try:
    assert len(transform_D4xD4)==96
    # save recomputation in case of multiple re-loading of the file
except NameError:
    print ("starting to compute special transformations transform_D4xD4.2")

    special_poss1 = [(x,1) for x in (Q(1),i,j,k)]
    special_poss2 = [(x,1/2) for x in (1-k, i+j,i-j,j+k,j-k,i+k,i-k,
                                            1+i,1-i,1+j,1-j,1+k)] # empirically 1-k and i+j are sufficient
    special_poss3 = [(x,1/4) for x in (1+i+j+k,1+i+j-k,1+i-j+k,1+i-j-k,
                                       1-i+j+k,1-i+j-k,1-i-j+k,1-i-j-k)]
    
    # units == list(a*b for a in (quat_1,quat_minus1) for b in (quat_1,quat_i,quat_j,quat_k))
    units = list(a*b for a in (quat_1,quat_minus1) for b in (quat_1,quat_i,quat_j,quat_k)) ## ADDED 2022
    elements_D4xD4 = [(FDR(l,r,star),FDR(l.to_exact_quaternion(),r.to_exact_quaternion(),star))
                      for l in units for r in units for star in (False, True)]
    exact_to_2D = {ex:e for e,ex in elements_D4xD4}
    transf = set()
    transform_D4xD4 = []
    #good_conj = []
    for special_poss in (special_poss1,special_poss2,special_poss3,):
        # mixed between different does not produce valid solutions  
        for l,scale1 in special_poss:
            for r,scale2 in special_poss:
                res = dict()
                #print(l,r,"scales",scale1,scale2)
                try:
                    for e,e_exact in elements_D4xD4:
                        ee = AxB_group.conjugate_element(e_exact,l,r,scale1,scale2)
                        res[e]=exact_to_2D[ee]
                    transf.add(tuple(res[e] for e,_ in elements_D4xD4))
                    transform_D4xD4.append((res,(l,r,scale1,scale2)))
                except KeyError:
                    pass

    # rearrangement: the most frequently used transformations come first
    transform_D4xD4 = (transform_D4xD4[16:20] + transform_D4xD4[36:40] +
                       transform_D4xD4[:16]   + transform_D4xD4[20:36] + transform_D4xD4[40:])
    
    print ("done with special transformations transform_D4xD4")

    # special conjugations for one SO(3) group (left or right) individually,
    # for chiral groups.
    # these replace the "special" conjugations ("S",i) as part of the quat_2D class
    elements_D4 = [(l,l.to_exact_quaternion()) for l in units]
    exact_to_2D = {ex:e for e,ex in elements_D4}
    transf1 = set()
    transform_D4 = []
    for l,scale in (special_poss1[:1] # identity first
                    + special_poss2 + special_poss1[1:] + special_poss3):
         
        res = {e: exact_to_2D[l.conjugate()*e_exact*l*scale] for e,e_exact in elements_D4}
        transf1.add(tuple(res[e] for e,_ in elements_D4))
        transform_D4.append((res,(l,scale)))

assert len(transf1)==len(transform_D4)==24
assert len(transf)==len(transform_D4xD4)==96
# confirmed by more elaborate computations        
## 96 elements. What is this group?

# Actually, almost all of this is never used. only 8 of them involving 1-k, 1+k, i+j,and i+j!

