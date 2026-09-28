# 1-fromFSPT

Exact cochain calculations for fermionic phases with one-form symmetry
in 2+1D through 4+1D, and finite-Abelian classification calculators for
3+1D and 4+1D.

The code accompanies the one-form solutions in the manuscript
*Lattice Model of Higher-Form/Group Fermionic Symmetry-Protected Topological
Phases*. It uses the analytic obstruction and stacking convention of
Ning, Ren, Wang, Qi, and Gu, *Stacking Group of 3+1D Interacting Fermionic
Symmetry-Protected Topological Phases* (manuscript revised September 26, 2026).

## Install and run

Python 3.12 is the reference environment. The dependencies are NumPy and
SymPy; their tested versions are pinned in `requirements.txt`.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt

python cochain_21d_verify.py --output results/cochain_21d_checks.json
python cochain_oneform_verify.py --output results/cochain_oneform_checks.json
python verify_classification_snf.py --output results/classification_checks.json
```

With `uv`, replace the environment setup with:

```bash
uv venv --python 3.12
uv pip install -r requirements.txt
. .venv/bin/activate
```

## Calculate a 3+1D classification

Write the symmetry as a product of cyclic groups with orders `N_i` and
specify the parity character by coefficients `epsilon_i` in `{0,1}`.
A coefficient must be zero on an odd-order factor.

```bash
python classify.py --orders 2 --twist 1
python classify.py --orders 4,6 --twist 1,0
python classify.py --orders 8 --twist 1 --absolute-bordism
```

The output reports the invariant factors and group order. An empty list of
invariant factors denotes the trivial group. The first two examples give
`Z_16` and `Z_4 x Z_48`; the last gives `Z_2`.

For a split lift, supply all-zero character coefficients. For the trivial
symmetry, use `--orders 1 --twist 0`.

## Which 3+1D classification is being computed?

The default is the paper's **relative one-form classification**. Its data
are `(n_2, n_3, nu_4)`: Majorana-chain, complex-fermion, and bosonic
decorations. There is no symmetry-dependent `n_1` occupation in this pure
one-form construction. For a nonzero parity character, the relative
classification contains a canonical universal subgroup `Z_16`.

More precisely, write `G_w = Omega_4^w(B^2 A)` for absolute twisted-spin
bordism, and let `i` include ordinary spin bordism at the trivial
background. The default calculation gives

```text
F(A,w) = Hom(G_w / image i, U(1)).
```

`--absolute-bordism` computes a separate mathematical object,

```text
D(A,w) = Hom(Tor G_w, U(1)).
```

For nonzero `w`, the map from `F` to `D` takes the quotient by the
canonical universal `Z_16`; for zero `w`, the two groups agree. Within
this mathematical quotient, every class is represented by the image of
the bosonic sector. This comparison does **not** establish an additional
symmetric local `p+ip` equivalence of the lattice model. In particular,
the trivial absolute comparison for twisted `Z_2` and `Z_4` does not
replace their relative one-form classification `Z_16`.

The original `--full` name overstated that physical interpretation.
It is retained as a deprecated alias for `--absolute-bordism` and prints
a clarification to standard error. Both spellings report the same
mathematical comparison in JSON. See [CHANGELOG.md](CHANGELOG.md).

The scope is **finite Abelian, unitary, pure one-form symmetry** with a
parity character. Continuous symmetries and nontrivial two-group
Postnikov data are outside this calculation.

## 3+1D mathematical input

Let `X = B²A` and `w` be the parity character. The bosonic subgroup is

```text
B(A,w) = H^4(X,U(1)) / image R_w,
R_w(x) = (x^2 + w x)/2.
```

For nonzero `w`, adjoining the universal Majorana root `M` gives

```text
Sq^1 w != 0:  4M = q_w,    q_w = [-P_2(w)/4],
Sq^1 w == 0:  2M = r_w,    r_w = [-P_4(lambda)/8],
```

where `lambda: A -> Z_4` lifts `w`. The class `r_w` is independent of that
lift after taking the bosonic quotient. For `w=0`, the answer is simply
`B(A,0)`.

For an even cyclic group with nonzero character, the relative one-form group is

| Cyclic order | Classification |
| --- | --- |
| `N = 2 mod 4` | `Z_(8N)` |
| `N = 4 mod 8` | `Z_(4N)` |
| `N = 0 mod 8` | `Z_(2N) x Z_2` |

The scripts evaluate the full mixed quadratic presentation for arbitrary
finite products. They do not infer product-group answers by multiplying
single-factor classifications.

## 4+1D classification and no-go theorem

The full decoration data are `(n_2,n_3,n_4,nu_5)`, where `n_2` is the
integer-valued `p+ip` occupation, `n_3` the Majorana occupation, and `n_4`
the complex-fermion occupation. For finite Abelian, unitary, pure
one-form symmetry, every class has a bosonic representative, for every
parity character. This does not mean that every response is trivial.

```bash
python classify_41d.py --orders 2 --twist 1
python classify_41d.py --orders 4 --twist 1
python classify_41d.py --orders 2,2 --twist 1,0
python verify_41d_classification.py --output results/classification_41d_checks.json
python cochain_41d_verify.py --output results/cochain_41d_checks.json
```

These examples return `Z_2`, the trivial group, and `Z_2 x Z_2`.
The calculation starts with all four layers: `H^2(B²A,Z)=0` removes the
`p+ip` occupation, while `H^1(B²A,Z)=0` removes residual integer
equivalences. The Majorana layer is obstructed, and every unobstructed
complex-fermion occupation is an incoming Majorana coboundary, with its
compensating phase included. The bosonic quotient is

```text
F_(4+1)(A,w) = H^5(B²A,U(1)) / < (Sq^2 v + w v)/2 : v in H^3(B²A,Z_2) >.
```

For `A = product_i Z_(N_i)`, its abstract group is the product of
`Z_gcd(N_i,N_j)` over `i<j`, with one extra `Z_2` precisely when
`Sq^1 w != 0`. In cyclic coordinates this condition means that some
active character coefficient has `N_i = 2 mod 4`. The extra response
is the bosonic action `(w cup Sq^1 w)/2`, equal to `(w_2 cup w_3)/2`
on the allowed twisted-spin backgrounds.

The calculator takes the Smith normal form of the full cohomological
quotient in the supplied cyclic coordinates. It does not assume that
raw mixed generators in those coordinates survive independently.
`verify_41d_classification.py` compares this quotient with the abstract
product formula, and separately checks the complex-layer kernel/image
identity by binary linear algebra. The cochain checker verifies the
compensating Majorana phase and the vanishing secondary identification.
The manuscript supplies the cohomology and bordism proofs; finite
arithmetic checks alone do not establish completeness.

In degree five the absolute bordism characters, the quotient by the
ordinary-spin image, and the cofiber-relative convention agree. There
is therefore no separate `--absolute-bordism` option for this calculator.
The existing 3+1D calculator and its conventions are unchanged.

## 2+1D no-go checks

For finite Abelian, unitary, pure one-form symmetry, the manuscript proves
that every symmetry-dependent 2+1D response is trivial, for either a split
or twisted fermionic lift. Purely gravitational phases are excluded, and
a choice of trivialization at zero background is not additional phase data.

The obstruction permits only the complex-fermion occupations `n_2=0` and,
for nonzero twist `w`, `n_2=w`. A contractible Kitaev-loop equivalence
identifies them. With additive phases, its cochain expression is

```text
(n_2, nuhat_3) ~ (n_2+w, nuhat_3 + (n_2 cup_1 w)/2).
```

The incoming operation is `H^0(X,Z_2) -> H^2(X,Z_2), 1 -> w` in the
binary Majorana coefficient row. Although `H^1(X,Z_2)=0` rules out a
Majorana occupation, it does not remove this equivalence. The microscopic
argument uses auxiliary vacuum pairs, triangle Kitaev loops, and local
reconnections with the one-form transport signs. Its separate
cohomological formulation uses the twisted differential `Sq^2 + w cup`.

`cochain_21d_verify.py` checks the phase identities with integer arithmetic:

- All 4,096 pairs of independent binary occupation and twist cocycles on
  a four-simplex for compatibility with the obstruction.
- The square of the equivalence on those pairs, verifying that its
  remaining phase is the coboundary `d tilde(w)/4`.
- All eight binary two-cocycles on a tetrahedron for the explicit
  trivialization of the candidate `n_2=w`.

All three residual counts vanish. These are cochain checks; the script
does not reconstruct the auxiliary-fermion circuit or prove completeness.
It imports the higher-cup operations from `cochain_oneform_verify.py`.

## 3+1D cochain checks

`cochain_oneform_verify.py` implements normalized interval-cut products
and signed integer higher cups. All phase comparisons use integer
numerators modulo 8 or 16.

The current four-word Adem representative is essential. With
`U = tilde(u)`, `P = U cup U + U cup_1 dU`, `b = u cup_1 u`, and
`y = b cup_1 u`, the chosen roots have top cochains

```text
B: -P/4,        C: -P/8,        M: 3P/16 + y/2.
```

Their stacking relations are `2M=C`, `2C=B`, and `4B=0`.
The script checks:

- All 1,024 binary five-simplices for the root equations and stacking
  primitives.
- All 1,048,576 `Z_4` five-simplices for the exact complex-occupation
  removal and its bosonic reduction.
- All 262,144 `Z_4 x Z_2` four-simplices for independence of the character
  lift modulo `R_w`.
- The restricted microscopic-to-analytic coordinate dictionary on 1,024
  obstruction inputs and 64 identical-pair inputs.

The earlier `P/16 + y/2` root is an intentional negative control: it fails
on 212 binary five-simplices in the current convention. This expected
count is separate from the zero residuals of all corrected formulas.

`cochain_microscopic_diagonal.json` contains the finite microscopic input
arrays and coordinate witnesses with provenance hashes. Those arrays are
fixed inputs from the zero-form calculation. The checker verifies their
coordinate dictionary; it does not regenerate the original projector
amplitudes.

## Classification checks

`verify_classification_snf.py` compares two constructions:

1. Smith reduction of the full quadratic and fermionic presentation in
   the original cyclic coordinates.
2. The closed product formula after primary decomposition and a change
   of basis adapted to the character.

The reference run checks 5,858 presentations, covering the relative
classification and the absolute-bordism comparison, all characters on up to three two-primary factors with
exponents 1 through 4, cyclic orders through 32, and a family of odd and
mixed-order products. Eight invalid inputs are also checked.
All comparisons pass. These arithmetic checks supplement the
cohomological proof; they do not establish its topological inputs.

Recorded results are in `results/`. The classification report includes
the script hash, Python and SymPy versions, timestamp, and coverage.
GitHub Actions runs the cochain and classification checkers in all three
dimensions and the 3+1D CLI interpretation tests. Run the latter locally with
`python -m unittest test_cli.py`.

## References

The mathematical inputs include:

- H. Cartan, [mod-two cohomology of Eilenberg-Mac Lane spaces](https://numdam.org/item/SHC_1954-1955__7_1_A10_0/)
  and [integral homology](https://numdam.org/item/SHC_1954-1955__7_1_A11_0/),
  Séminaire Henri Cartan 7 (1954–1955), Exposés 10 and 11.
- A. Kapustin and R. Thorngren,
  [Topological Field Theory on a Lattice, Discrete Theta-Angles and Confinement](https://arxiv.org/abs/1308.2926),
  for universal quadratic groups and Pontryagin squares.
- Q.-R. Wang and Z.-C. Gu,
  [Construction and Classification of Symmetry-Protected Topological Phases in Interacting Fermion Systems](https://arxiv.org/abs/1811.00536),
  for the auxiliary Kitaev-loop equivalence in the zero-form construction.
  The manuscript checks the transport signs for its one-form adaptation.
- T. Johnson-Freyd,
  [(3+1)D Topological Orders with Only a Z2-Charged Particle](https://arxiv.org/abs/2011.11165),
  and R. Kobayashi, A. Prem, and M. Yu,
  [Sixteen-Fold Way for Fermionic Topological Orders](https://arxiv.org/abs/2606.28682),
  for twisted supercohomology and the relative `Z_16` convention.
  Johnson-Freyd also gives the explicitly one-form differential used in
  the 2+1D no-go argument.
- D. S. Freed and M. J. Hopkins,
  [Reflection Positivity and Invertible Topological Phases](https://arxiv.org/abs/1604.06527),
  for the Anderson-dual universal coefficient sequence used in the
  absolute-bordism comparison. This does not by itself establish an extra
  higher-form lattice equivalence.
- A. Debray and M. Yu,
  [What Bordism-Theoretic Anomaly Cancellation Can Do for U](https://arxiv.org/abs/2210.04911),
  Corollary 4.44 and Proposition 4.45, for the Wu-manifold characteristic
  number used to exclude the remaining 4+1D differential.
