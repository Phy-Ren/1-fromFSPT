# Changelog

## Unreleased — classification scope clarification

The initial release called the absolute-bordism comparison a "full"
classification obtained by allowing integer `p+ip` redefinitions. That
wording asserted a physical equivalence of the one-form lattice model
that had not been established. The mathematical comparison and the
calculated groups remain unchanged.

- The default is explicitly the paper's relative one-form classification
  with data `(n_2, n_3, nu_4)`. Its universal `Z_16` subgroup is retained.
- `--absolute-bordism` computes characters of the torsion of absolute
  twisted-spin bordism. It is presented as a mathematical comparison,
  without claiming an additional allowed lattice move.
- `--full` remains a deprecated CLI alias and emits an explanatory warning.
- JSON convention and scope fields, documentation, verifier provenance,
  and CI examples now use these meanings. Python helper arguments use
  `absolute_bordism` in place of the former `full` keyword.
- The integer presentations, cochain formulas, and classification
  invariant factors have not changed.
