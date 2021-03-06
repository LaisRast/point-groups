# 4-Dimensional Point Groups

This repository contains scripts written in [Sage](https://www.sagemath.org/)
for the article [Towards a Geometric Understanding of the 4-Dimensional
Point Groups](http://arxiv.org/abs/2205.04965) by Laith Rastanawi and Günter Rote.

## Overview

* `logs/`: Contains:
  * Outputs of the scripts `classify_AxA.sage` and `classify_DxD.sage`:
    `TxT.log`, `OxO.log`, `IxI.log` and `DxD_{orderbound}.log`.
  * Cross-reference between CS names and
    our names for the toroidal groups:
    `toroidal_chiral_cross_reference.log`
    and
    `toroidal_achiral_cross_reference.log`.

* `precomputed/`:
  Contains precomputed data such as:
  * Groups catalogs.
  * Order catalogs.

* `main_functions.sage`:
  Contains definitions, classes and functions that are used in the
  other scripts. Most notable ones:
  * `Q([x1, x2, x3, x4])`:
    Represents the quaternion $x_1 + x_2 i + x_3 j + x_4 k$.
  
  * `quat_2D(alpha, p)`:
    Represents the quaternion $e^{\alpha \pi} j^p$,
    an element of the (binary) dihedral group $2D_2n$.
  
  * `FDR(l, r, star)`:
    Represents the 4-dimensional orthogonal transformation
    $[l, r]$ if `star == False`, and $\star[l,r]$ otherwise.
  
  * `Q_group(gens)`:
    Represents the quaternion group generated by `gens`.
  
  * `AxB_group(gens)`:
    Represents the 4-dimensional point group generated by `gens`.
  
  * `classify_group(gr)`:
    Returns the fingerprint of the group `gr`.

  * `Hurley_pattern(fingerprint)`:
    Converts the fingerprint `fingerprint` of a point group into Hurley pattern.
  
  * `mk_group(gens)`:
    Returns the group generated by `gens`.
  
  * `toroidal_group(*parameters)`:
    Returns the toroidal group defined by the parameters `parameters`.

  * Functions to generate catalogs of point groups:
    * `generate_polyhedral_and_axial` 
    * `generate_chiral_toroidal`
    * `generate_achiral_toroidal`
    * `generate_tubular`

* `generate_toroidal_CS.sage`:
  Contains functions to generate catalogs for the
  toroidal groups using CS classification.
 
* `generate_catalog.sage`:
  Generates groups catalogs and stores them in `precomputed/` as json files.
  The keys are fingerprints of the groups.
  The values are the names of the groups.
  
* `counting_groups.sage`:
  Generates order catalogs and stores them in `precomputed/` as json files.
  The keys are orders of the groups.
  The values are the names of the groups.
  
* `generate_oeis_sequences.sage`:
  Generates the following two sequences to be used in OEIS:
  * number of 4-dimensional point groups of order n, and
  * number of 4-dimensional chiral point groups of order n.
  
  The precomputed catalogs are generated using the main 
  function in `counting_groups.sage` with `orderbound = 10000`.
 
* `classify_AxA.sage`:
  Classifies all subgroups of $\pm[A\times A]$, 
  and their achiral extensions.
  Outputs the results in `logs/`.

* `classify_DxD.sage`:
  Classifies all chiral toroidal groups up to order `orderbound//2`
  and the achiral toroidal groups up to order `orderbound`.
  Outputs the results in `logs/`.

* `conjugacy_Dn.sage`:
  Generates dictionaries `transform_D4` and `transform_D4xD4` to compute
  conjugations for elements of $2D_{2n}$ by non-elements of $2D_{2n}$.

* `table_toroidal_duplications.sage`:
  Checks if we have the correct conjugations for toroidal groups and
  generates a LaTeX table `tables/toroidal_duplications.tex` of them.

* `table_polyhedral_groups.sage`:
  This script has two parts:
  * Checks that our generators and CS generators for the polyhedral and the
    axial groups generate geometrically equivalent (but not necessarily equal)
    groups.
    The list of generators is given in `polyhedral_and_axial_lists.sage`.
  * Generates 3 LaTeX tables:
      * `tab_polyhedral_and_axial_check.tex`:
        For internal use, with the full information.
      * `tab_polyhedral_and_axial_appendix.tex`:
        All polyhedral and axial groups, with generators.
      * `tab_polyhedral.tex`:
        An overview of the polyhedral groups.

* `table_axial_groups.sage`:
  This script generates 4 LaTeX tables:
  * `tab_axial_combined.tex`:
    For internal use, with the full information.
  * `tab_axial_full.tex`:
    All axial groups with cross-references.
  * `tab_axial_pap.tex`:
    Pyramidal and prismatic.
  * `tab_axial_hybrid.tex`:
    Hybrid groups with methods.
    
  This files loads `table_polyhedral_groups.sage` to generate a catalog.

* `polyhedral_and_axial_lists.sage`:
  Contains two lists of polyhedral and axial groups:
  * `polyhedral_list_Cox` whose entries consist of 
  [DuVal number,  DuVal name, CS generators, our genreators, Coxeter name, method]
  * `axial_list` whose entries consist of
  [DuVal number,  DuVal name, CS generators, our genreators]

* `axial_data.sage`:
  Contains data used in `table_axial_groups.sage`.
  The Coxeter names in this files are
  also used in `table_polyhedral_groups.sage`.

* `BBNWZ_cat.sage`:
  Contains the dictionary `BBNWZ_cat`
  which is used to convert to numbering of
  crystallographic point groups.
  
* `generate_tubes.sage`:
  Generates decompositions of polar orbit polytopes of cyclic-type tubical
  groups $G$ when the image of the starting point is a rotation center of $G^h$.
  Outputs are json files which are stored in `data/`.
  These files are used in the online gallery.
