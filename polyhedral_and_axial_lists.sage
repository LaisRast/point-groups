#!/usr/bin/env python 
# coding: utf-8

## Contains two lists of polyhedral and axial groups:
## * `polyhedral_list_Cox` whose entries consist of 
## [DuVal number,  DuVal name, CS generators, our genreators, Coxeter name, method]
## * `axial_list` whose entries consist of
## [DuVal number,  DuVal name, CS generators, our genreators]
## DuVal number is only for information and is not checked.
## All diploid CS generators should contain "-".

polyhedral_list_Cox = [
    (50, r"(I/I;I/I)^*",
     [[i_I,1],[w,1],[1,i_I],[1,w],";","-","*"],
     [[i_I,w],[w,i_I], "*"], "[3,3,5]", ""),	
    (30, r"(I/I;I/I)",
     [[i_I,1],[w,1],[1,i_I],[1,w],";","-"],
     [[i_I,w],[w,i_I]], "[3,3,5]^+", "chiral part"),	
    (29, "(I/I;O/O)",
     [[i_I,1],[w,1],[1,i_O],[1,w],";","-"],
     [[i_I,w],[w,i_O]], r"[[3,3,5]^+_{\frac15L}]", r"inscribed polar \& swap"), #$ '+-[IxO]
    (29, "(O/O;I/I)",
     [[i_O,1],[w,1],[1,i_I],[1,w],";","-"],
     [[i_O,w],[w,i_I]], r"[[3,3,5]^+_{\frac15R}]", r"inscribed polar \& swap"), #$ '+-[OxI]
    (24, "(I/I;T/T)",
     [[i_I,1],[w,1],[1,i],[1,w],";","-"],
     [[i_I,w],[w,i]], r"[3,3,5]^+_{\frac15L}", "inscribed polar"), #$ '+-[IxT]
    (24, "(T/T;I/I)",
     [[i,1],[w,1],[1,i_I],[1,w],";","-"],
     [[i,w],[w,i_I]], r"[3,3,5]^+_{\frac15R}", "inscribed polar"), #$ '+-[TxI]
    (48, "(O/O;O/O)^*", #$&   $[...
     [[i_O,1],[w,1],[1,i_O],[1,w],";","-","*"],
     [["*",1,w],[w,i_O]], "[[3,4,3]]", ""), #$ '+-[OxO].2'
    (25, "(O/O;O/O)", #$&   $[...
     [[i_O,1],[w,1],[1,i_O],[1,w],";","-"],
     [[i_O,w],[w,i_O]], "[[3,4,3]]^+", "chiral part"), #$ '+-[OxO]
    (45, "(O/T;O/T)^*",
     [[i,1],[w,1],[1,i],[1,w],";","-",[i_O,i_O],"*"],
     [["*",w,1],[i_O,i_O]], "[3,4,3]", "nonswapping"),
    (46, "(O/T;O/T)^*_-",
     [[i,1],[w,1],[1,i],[1,w],";","-",[i_O,i_O],["*",1,i_O]],
     [[w,1],["*",1,i_O]], "[[3,4,3]^+]", "swap with mirror"),
    (28, "(O/T;O/T)", #$&
     [[i,1],[w,1],[1,i],[1,w],";","-",[i_O,i_O]],
     [[w,1],[1,w],[i_O,i_O]], "[3,4,3]^+", r"chiral \& nonswapping"), #$ '+-1/2[OxO]
    (43, "(T/T;T/T)^*", #$&   $[3,4,3^+]$
     [[i,1],[w,1],[1,i],[1,w],";","-","*"],
     #[[i,w],["*",i_O,i_O]]), 
     [[i,w],["*",w,i]], "[3,4,3^+]", "edge orientation"), #$ %1152 '+-[TxT].2'
    (23, "(O/O;T/T)", #$&   $[...
     [[i_O,1],[w,1],[1,i],[1,w],";","-"],
     [[i_O,w],[w,i]], "[[^+3,4,3^+]]_L", "diagonal marking"), #$ '+-[OxT]
    (23, "(T/T;O/O)", #$&   $[...
     [[i,1],[w,1],[1,i_O],[1,w],";","-"],
     [[i,w],[w,i_O]], "[[^+3,4,3^+]]_R", "diagonal marking"), #$ '+-[TxO]
    (20, "(T/T;T/T)",
     [[i,1],[w,1],[1,i],[1,w],";","-"],
     [[i,w],[w,i]], "[^+3,4,3^+]", "2 dual edge orientations"), #$ '+-[TxT]
    (47, "(O/V;O/V)^*",
     [[i,1],[j,1],[1,i],[1,j],";","-",[w,w],[i_O,i_O],"*"],
     #[["*",i,1],[w,w],[i_O,i_O]], "", ""),
     [["*",i*w,w],[i_O,i_O]], "[3,3,4]", ""),
    (27, "(O/V;O/V)", #$&#$ '+-1/6[OxO]  
     [[i,1],[j,1],[1,i],[1,j],";","-",[w,w],[i_O,i_O]],
     [[i,j],[w,w],[i_O,i_O]], "[3,3,4]^+", "chiral part"),
    (41, "(T/V;T/V)^*", # 41/42 numbering according to Goursat; the description in Du Val is confused.
     [[i,1],[j,1],[1,i],[1,j],";","-",[w,w],"*"],
     [["*",i,1],[w,w]], "[^+3,3,4]", "even permutations"), 
    (42, "(T/V;T/V)^*_-",
     [[i,1],[j,1],[1,i],[1,j],";","-",[w,w_bar],"*"],
     [["*",i,1],[w,w_bar]], "[3,3,4^+]", "2-coloring"),
    (22, "(T/V;T/V)", #$\pm \frac13[T\times T]$
     [[i,1],[j,1],[1,i],[1,j],";",[w,w],"-"],
     [[i,1],[1,i],[w,w]], "[^+3,3,4^+]", r"2-coloring \& chiral"),
    (51, r"(I^\dag/C_2;I/C_2)^{\dag* }",
     [";",[w,w],[i_I,i_I2], "-","*"],
     [[w,-w],["*",i_I,i_I2]], "[[3,3,3]]", ""),
    (51, r"(I^\dag/C_2;I/C_2)^{\dag*}",
     [";",[w,w],[i_I,i_I2], "-","*"],
     [[w,-w],["*", i_I*i_O*i, i_I_dag*i_O*i]], "[[3,3,3]]", ""),
    (32, r"(I^\dag/C_2;I/C_2)^{\dag}",
     [";",[w,w],[i_I,i_I2],"-"],
     [[w,w],[i_I,-i_I2]], "[[3,3,3]]^+", "chiral part"),
    (32, r"(I^\dag/C_2;I/C_2)^{\dag}",
     [";",[w,w],[i_I,i_I2],"-"],
     [[w,w],[i_I,-i_I_dag]], "[[3,3,3]]^+", "chiral part"),
    ("51'", r"(I^\dag/C_1;I/C_1)^{\dag*}", # "+1/60[IxI'].2_1"
     [";",[w,w],[i_I,i_I2], "-*"],
     [[w,w],[i_I,i_I2], "-*"], "[3,3,3]", "nonswapping"),
    ("51'", r"(I^\dag/C_1;I/C_1)^{\dag*}", # "+1/60[IxI'].2_1"
     [";",[w,w],[i_I,i_I2], "-*"],
     [[w,w],["*", i_I*i_O*i, i_I_dag*i_O*i]], "[3,3,3]", "nonswapping"),
    ("51'", r"(I^\dag/C_1;I/C_1)^{\dag*}_ -",
     [";",[w,w],[i_I,i_I2], "*"],
     [[w,w],["*", i_I,i_I2]], r"[[3,3,3]^+]", "swap with mirror"), ## ALTERNATIVE: "[3,3,3]^\circ"
    ("51'", r"(I^\dag/C_1;I/C_1)^{\dag*}_-", # "+1/60[IxI'].2_3"
     [";",[w,w],[i_I,i_I2], "*"],
     [[w,w],["*", i_I*i_O*i, -i_I_dag*i_O*i]], r"[[3,3,3]^+]", "swap with mirror"), ## ALTERNATIVE: "[3,3,3]^\circ"
    ("32'", r"(I^\dag/C_1;I/C_1)^{\dag}",
     [";",[w,w],[i_I,i_I2]],
     [[w,w],[i_I,i_I2]], "[3,3,3]^+", r"chiral \& nonswapping"),
    ("32'", r"(I^\dag/C_1;I/C_1)^{\dag}",
     [";",[w,w],[i_I,i_I2]],
     [[w,w],[i_I,i_I_dag]], "[3,3,3]^+", r"chiral \& nonswapping"),
]

axial_list = [
    # pyramidal and prismatic
    ("49'", "(I/C_1;I/C_1)^*", # & $[3,5]$ & $120$  .2_3
     [";",[w,w],[i_I,i_I],"*"],
     [[w,w],["*",i_I,i_I]]),
    (49, "(I/C_2;I/C_2)^*", # $2.[3,5]=[2,3,5]$ & $240$\\ &$\pm\frac1{60}[I\times I]\cdot 2$
     [";",[w,w],[i_I,i_I],"-","*"],
     [[w,-w],["*",i_I,-i_I]]),
    # $+ I$ &$+\frac1{60}[I\times I]$
    ("31'", "(I/C_1;I/C_1)", #$& $[3,5]^+$ & $60$
     [";",[w,w],[i_I,i_I]],
     [[w,w],[i_I,i_I]],),
    ("49'", "(I/C_1;I/C_1)^*_-",# $[3,5]^\circ$ & $120$\\ &$+\frac1{60}[I\times I]\cdot 2_1$
     [";",[w,w],[i_I,i_I],"-*"],
     [[w,w],["*",i_I,-i_I]]),
    # $\pm O$ &$+\frac1{24}[O\times O]\cdot 2_3$
    ("44'", r"(O/C_1;O/C_1)^{*\prime}",# $& $[3,4]$ & $48$ &
     [";",[w,w],[i_O,i_O],"*"],
     [[w,w],[i_O,i_O],"*"],),
    (44, "(O/C_2;O/C_2)^*", # $2.[3,4]=[2,3,4]$ & $96$$\pm\frac1{24}[O\times O]\cdot 2$
     [";",[w,w],[i_O,i_O],"-","*"],
     [[w,-w],[i_O,i_O],"-*"],),
    # $+O$ &$+\frac1{24}[O\times O]$
    ("26'", "(O/C_1;O/C_1)'", # $[3,4]^+$ & $24$ 
     [";",[w,w],[i_O,i_O]],
     [[w,w],[i_O,i_O]],),
    ("44'", r"(O/C_1;O/C_1)^{*\prime}_-", # $[3,4]^\circ$ & $48$\\$+\frac1{24}[O\times O]\cdot 2_1$
     [";",[w,w],[i_O,i_O],"-*"],
     [[w,w],[i_O,i_O],"-*"],),
    # $TO$ &$+\frac1{12}[T\times \overline T]\cdot 2_1$
    ("40'", "(T/C_1;T/C_1)^*  ", # $[3,3]$ & $24$ &
     [";",[w,w_bar],[i,-i],"-*"],
     [[w,w_bar],["*",i,i]]),
    ("40'", "(T/C_1;T/C_1)^* ", # $[3,3]$ & $24$ &
     [";",[w,w_bar],[i,-i],"-*"],
     [[w,w],["*",i_O,i_O]]),
    ("44''", r"(O/C_1;O/C_1)^{*\prime\prime}_-", #$& $[2,3,3]$ & $48$\\$+\frac1{24}[O\times \overline O]\cdot 2_1$
     [";",[w,w],[i_O,-i_O],"-*"],
     [[w,w],[i_O,-i_O],"-*"],),
    # $\pm T$ &$+\frac1{12}[T\times T]\cdot 2_3$&
    ("39'", "(T/C_1;T/C_1)^*_c", # & $[^+3,4]$ & $24$  
     [";",[w,w],[i,i],"*"],
     [[w,w],["*",i,i]],),
    (39, "(T/C_2;T/C_2)^*_c", # $& $2.[^+3,4]$ & $48$\\&$\pm\frac1{12}[T\times T]\cdot 2$
     [";",[w,w],[i,i],"-","*"],
     [[w,-w],["*",i,-i]],),
    #  $+ T$ &$+\frac1{12}[T\times T]$
    ("21'", "(T/C_1;T/C_1)", #$& $[3,3]^+$ & $12$ 
     [";",[w,w],[i,i]],
     [[w,w],[i,i]],),
    ("39'", "(T/C_1;T/C_1)^*_{c-}", # & $[^+3,4]^\circ$ & $24$\&$+\frac1{12}[T\times T]\cdot 2_1$&
     [";",[w,w],[i,i],"-*"],
     [[w,w],["*",i,-i]],),
    # hybrid axial
    # $+I$ in $\pm I$ &$\pm\frac1{60}[I\times I]$&
    (31, "(I/C_2;I/C_2)", #   &$2.[3,5]^+$ & $120$ & center, chirality\\
     [";",[w,w],[i_I,i_I],"-"],
     [[w,w],[i_I,-i_I]],),
    #  $\pm T$ in $\pm O$ &$+\frac1{24}[O\times \overline O]\cdot 2_3$
    ("44''", r"(O/C_1;O/C_1)^{*\prime\prime}", # $[2,3,3]^\circ$ & $48$&edge orientation\\
     [";",[w,w],[i_O,-i_O],"*"],
     [[w,w],[i_O,-i_O],"*"],),
    #$+O$ in $\pm O$ &$\pm\frac1{24}[O\times O]$&
    (26, "(O/C_2;O/C_2)", #$& $2.[3,4]^+$ & $48$& center, chirality\\
     [";",[w,w],[i_O,i_O],"-"],
     [[w,-w],[i_O,i_O]],),
    #$ TO$ in $\pm O$ &$\pm\frac1{12}[T\times\overline T]\cdot 2$&
    (40, "(T/C_2;T/C_2)^ *", #$& $2.[3,3]$ & $48$&center, alternation\\
     [";",[w,w_bar],[i,-i],"-","*"],
     [[w,-w_bar],["*",i,-i]]),
    (40, "(T/C_2;T/C_2)^*", #$& $2.[3,3]$ & $48$&center, alternation\\
     [";",[w,w_bar],[i,-i],"-","*"],
     [[w,-w],["*",i_O,-i_O]]),
    # $+ T$ in $\pm T$ &$\pm\frac1{12}[T\times T]$
    (21, "(T/C_2;T/C_2)", # $2.[3,3]^+$ & $24$&center, chirality\\
     [";",[w,w],[i,i],"-"],
     [[w,-w],[i,i]]),
    # $+ T$ in $+O$ &$+\frac1{12}[T\times\overline T]\cdot 2_3$&
    ("40'", "(T/C_1;T/C_1)^*_ -", #& $[3,3]^\circ$ & $24$&alternation\\
     [";",[w,w_bar],[i,-i],"*"],
     [[w,w_bar],["*",i,-i]]),
    ("40'", "(T/C_1;T/C_1)^*_-", #& $[3,3]^\circ$ & $24$&alternation\\
     [";",[w,w_bar],[i,-i],"*"],
     [[w,w],["*",i_O,-i_O]]),
    # $+ T$ in $ TO$ &$+\frac1{24}[O\times\overline O]$&
    ("26''", "(O/C_1;O/C_1)''", # $[2,3,3]^+$ & $24$& chirality\\\hline
     [";",[w,w],[i_O,-i_O]],
     [[w,w],[i_O,-i_O]],),
]
