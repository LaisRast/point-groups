#!/usr/bin/env python 
# coding: utf-8

## Contains data used in `table_axial_groups.sage`.
## The Coxeter names in this files are
## also used in `table_polyhedral_groups.sage`.

poly3 = (
    "+-I",
    "+I",
    "+-O",
    "+O",
    "TO",
    "+-T",
    "+T",
)

pyramidal = {
    "+-I": ("+1/60[IxI].2_3", '[3,5]'),
    "+I":  ("+1/60[IxI]",     '[3,5]^+'),
    "+-O": ("+1/24[OxO].2_3", '[3,4]'),
    "+O":  ("+1/24[OxO]",     '[3,4]^+'),
    "TO":  ("+1/12[TxT'].2_1", '[3,3]'),
    "+-T": ("+1/12[TxT].2_3",  '[^+3,4]'),
    "+T":  ("+1/12[TxT]",     '[3,3]^+')
}

prismatic = {
    "+-I": ("+-1/60[IxI].2",  r"2.[3,5]"),
    "+I": ("+1/60[IxI].2_1",  r"[3,5]^\circ"),
    "+-O": ("+-1/24[OxO].2",  r"2.[3,4]"),
    "+O": ("+1/24[OxO].2_1",  r"[3,4]^\circ"),
    "TO": ("+1/24[OxO'].2_1", r"[2,3,3]"),
    "+-T": ("+-1/12[TxT].2",  r"2.[^+3,4]"),
    "+T": ("+1/12[TxT].2_1",  r"[^+3,4]^\circ"),
}

hybrid = {
    ("+I", "+-I"):  ("+-1/60[IxI]",      r"2.[3,5]^+", "center, chirality"),
    ("+-T", "+-O"): ("+1/24[OxO'].2_3",  r"[2,3,3]^\circ", "edge orientation"),
    ("+O", "+-O"): ("+-1/24[OxO]",      r"2.[3,4]^+", "center, chirality"),
    ("TO", "+-O"): ("+-1/12[TxT'].2",   r"2.[3,3]", "center, alternation"),
    ("+T", "+-T"):  ("+-1/12[TxT]",      r"2.[3,3]^+", "center, chirality"),
    ("+T", "+O"):   ("+1/12[TxT'].2_3",  r"[3,3]^\circ", "alternation"),
    ("+T", "TO"):   ("+1/24[OxO']",      r"[2,3,3]^+", "chirality"),
}

poly3groups = {
    # string: (latex, orbitope, IT)
    # IT uses slanted "m"
    "+-I": (r"\pm I", "*532", "53m", 120),
    "+I":  ("+I", "532", "532", 60),
        # }$ &$+\frac1{60}[I\times I]$& $[3,5]^+$ & $60$
    "+-O": (r"\pm O", "*432", "m3m", 48),
        # }$ &$+\frac1{24}[O\times O]\cdot 2_3$& $[3,4]$ & $48$
    "+O": ("+O", "432", "432", 24),
        # }$ &$+\frac1{24}[O\times O]$& $[3,4]^+$ & $24$
    "TO": ("TO", "*332", r"\bar43m", 24),
        #   }$ &$+\frac1{12}[T\times \overline T]\cdot 2_1$& $[3,3]$ & $24$
    "+-T": (r"\pm T", "3{*}2", "m3", 24),
        # }$ &$+\frac1{12}[T\times T]\cdot 2_3$& $[^+3,4]$ & $24$
    "+T": ("+T", "332", "23", 12),
        #   }$ &$+\frac1{12}[T\times T]$& $[3,3]^+$ & $12$
}

poly3_generators = {
    "+-I":    [[w, w], [i_I, i_I], "*"],
    "+I":     [[w, w], [i_I, i_I]],
    "+-O":    [[w, w], [i_O, i_O], "*"],
    "+O":     [[w, w], [i_O, i_O]],
    "TO":     [[w, w], [i, i], ["*", i_O, i_O]],
    "+-T":    [[w, w], [i, i], "*"],
    "+T":     [[w, w], [i, i]],
}
