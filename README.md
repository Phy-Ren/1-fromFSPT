# 1-fromFSPT

Exact cochain calculations and finite-Abelian classification for 3+1D
fermionic phases with one-form symmetry.

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

python cochain_oneform_verify.py --output results/cochain_oneform_checks.json
python verify_classification_snf.py --output results/classification_checks.json
```

With `uv`, replace the environment setup with:

```bash
uv venv --python 3.12
uv pip install -r requirements.txt
. .venv/bin/activate
```

## Calculate a classification

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

## Which classification is being computed?

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

## Mathematical input

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

## Exact cochain checks

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
GitHub Actions runs the same two exact checkers and the CLI interpretation
tests. Run those tests locally with `python -m unittest test_cli.py`.

## References

The mathematical inputs include:

- H. Cartan, [mod-two cohomology of Eilenberg-Mac Lane spaces](https://numdam.org/item/SHC_1954-1955__7_1_A10_0/)
  and [integral homology](https://numdam.org/item/SHC_1954-1955__7_1_A11_0/),
  Séminaire Henri Cartan 7 (1954–1955), Exposés 10 and 11.
- A. Kapustin and R. Thorngren,
  [Topological Field Theory on a Lattice, Discrete Theta-Angles and Confinement](https://arxiv.org/abs/1308.2926),
  for universal quadratic groups and Pontryagin squares.
- T. Johnson-Freyd,
  [(3+1)D Topological Orders with Only a Z2-Charged Particle](https://arxiv.org/abs/2011.11165),
  and R. Kobayashi, A. Prem, and M. Yu,
  [Sixteen-Fold Way for Fermionic Topological Orders](https://arxiv.org/abs/2606.28682),
  for the twisted supercohomology and relative `Z_16` convention.
- D. S. Freed and M. J. Hopkins,
  [Reflection Positivity and Invertible Topological Phases](https://arxiv.org/abs/1604.06527),
  for the Anderson-dual universal coefficient sequence used in the
  absolute-bordism comparison. This does not by itself establish an extra
  higher-form lattice equivalence.
