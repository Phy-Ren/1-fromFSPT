# Changelog

## Unreleased — 4+1D finite one-form classification

- Added `classify_41d.py` for the full finite-Abelian bosonic quotient
  after all four decoration layers are included. It reports the absence
  of intrinsically fermionic classes without equating this with a trivial
  response group.
- Added independent Smith-form/product-formula comparisons and binary
  kernel/image checks: 6,254 presentations, 63 binary twists, and eight
  invalid inputs pass.
- Added exact cochain checks of the incoming Majorana coboundary, its
  compensating phase, complex stacking, and the secondary identification.
  All residuals vanish. Recorded outputs are in `results/`.
- Added the new checks to CI. Existing 2+1D and 3+1D algorithms and
  classification conventions are unchanged.

## Unreleased — 2+1D one-form checks

- Added `cochain_21d_verify.py` and its recorded results for the exact
  Kitaev-loop cochain identities supporting the finite-Abelian 2+1D
  no-go theorem. All 4,096 occupation/twist pairs and eight tetrahedral
  backgrounds pass.
- Documented the theorem's response convention, the separate physical
  equivalence argument, and the numerical checker's scope; added its
  run command to CI.
- The 3+1D algorithms, classification conventions, and numerical results
  are unchanged.

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
