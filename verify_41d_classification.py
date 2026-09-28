#!/usr/bin/env python3
"""Verify the 4+1D one-form quotient and absence of fermionic layers.

The integer check compares the full H^6 presentation in the original
cyclic coordinates against the independent product formula. The binary
check verifies exactness at quadratic polynomials of D_w = D + w,
D(u_i)=u_i^2, using F_2 linear algebra. The identification of these
polynomials with the obstruction modulo integral reductions is a
mathematical input; these tests do not replace the cohomology proof.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement, product
import json
from math import gcd
from pathlib import Path
import platform

import sympy
from classify_41d import presentation, validate_inputs


def combine_cyclic_factors(orders):
    """Align prime powers into invariant factors without Smith reduction."""
    powers = {}
    for n in orders:
        for prime, exponent in sympy.factorint(n).items():
            powers.setdefault(int(prime), []).append(int(exponent))
    length = max((len(values) for values in powers.values()), default=0)
    result = [1] * length
    for prime, exponents in powers.items():
        aligned = [0] * (length - len(exponents)) + sorted(exponents)
        for k, exponent in enumerate(aligned):
            result[k] *= prime ** exponent
    return tuple(result)


def closed_formula(orders, twist):
    validate_inputs(orders, twist)
    factors = [gcd(a, b) for a, b in combinations(orders, 2)]
    if any(w and n % 4 == 2 for n, w in zip(orders, twist)):
        factors.append(2)
    return combine_cyclic_factors(factors)


def gf2_rank(columns):
    """Rank a binary matrix whose columns are encoded as integer bit masks."""
    pivots = {}
    for column in columns:
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
    return len(pivots)


def polynomial_differential(monomial, twist, target_index):
    """Compute (D+w) on a monomial; a tuple records its variable indices."""
    result = 0
    for i, multiplicity in Counter(monomial).items():
        if multiplicity % 2:
            target = tuple(sorted(monomial + (i,)))
            result ^= 1 << target_index[target]
    for i, w in enumerate(twist):
        if w:
            target = tuple(sorted(monomial + (i,)))
            result ^= 1 << target_index[target]
    return result


def check_complex_layer(rank, twist):
    linear = list(combinations_with_replacement(range(rank), 1))
    quadratic = list(combinations_with_replacement(range(rank), 2))
    cubic = list(combinations_with_replacement(range(rank), 3))
    quadratic_index = {m: i for i, m in enumerate(quadratic)}
    cubic_index = {m: i for i, m in enumerate(cubic)}
    incoming = [polynomial_differential(m, twist, quadratic_index) for m in linear]
    obstruction = [polynomial_differential(m, twist, cubic_index) for m in quadratic]
    for column in incoming:
        composite = 0
        for j, image in enumerate(obstruction):
            if column & (1 << j):
                composite ^= image
        if composite:
            raise AssertionError({"rank": rank, "twist": twist, "failure": "D_w squared"})
    incoming_rank = gf2_rank(incoming)
    obstruction_rank = gf2_rank(obstruction)
    if incoming_rank + obstruction_rank != len(quadratic):
        raise AssertionError({"rank": rank, "twist": twist, "failure": "nonzero complex quotient"})
    if incoming_rank != rank - int(any(twist)):
        raise AssertionError({"rank": rank, "twist": twist, "failure": "wrong incoming kernel"})


def check_presentation(orders, twist):
    actual = presentation(orders, twist)
    expected = closed_formula(orders, twist)
    if actual != expected:
        raise AssertionError({"orders": orders, "twist": twist, "actual": actual, "expected": expected})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    counts = {"cyclic_presentations": 0, "product_presentations": 0,
              "trivial_presentations": 0, "binary_complex_layer_checks": 0,
              "invalid_inputs_rejected": 0}
    for n in range(1, 65):
        for w in range(1 if n % 2 else 2):
            check_presentation([n], [w])
            counts["cyclic_presentations"] += 1
    sample_orders = (1, 2, 3, 4, 6, 8, 9, 10, 12, 15, 16)
    for length in (2, 3):
        for orders in product(sample_orders, repeat=length):
            choices = [(0,) if n % 2 else (0, 1) for n in orders]
            for twist in product(*choices):
                check_presentation(orders, twist)
                counts["product_presentations"] += 1
    check_presentation([], [])
    check_presentation([1, 1, 1, 1], [0, 0, 0, 0])
    counts["trivial_presentations"] = 2
    for rank in range(6):
        for twist in product((0, 1), repeat=rank):
            check_complex_layer(rank, twist)
            counts["binary_complex_layer_checks"] += 1
    invalid = [([0], [0]), ([-2], [0]), ([2], []), ([3], [1]),
               ([2], [2]), ([True], [0]), ([2], [False]), ([2.0], [0])]
    for orders, twist in invalid:
        try:
            presentation(orders, twist)
        except ValueError:
            counts["invalid_inputs_rejected"] += 1
        else:
            raise AssertionError({"invalid_input_accepted": [orders, twist]})
    directory = Path(__file__).resolve().parent
    result = {
        "status": "passed", "failures": 0,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(), "sympy_version": sympy.__version__,
        "source_sha256": {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                          for name in ("classify_41d.py", "verify_41d_classification.py")},
        "counts": counts,
        "coverage": {
            "cyclic": "All valid twists for orders 1 through 64.",
            "products": "All ordered two- and three-factor tuples with orders 1,2,3,4,6,8,9,10,12,15,16 and all valid twists.",
            "complex_layer": "Every binary twist for zero through five degree-two generators; verifies D_w squared=0 and ker obstruction=image incoming.",
            "limitation": "Arithmetic checks of the quotient and polynomial exactness; topological inputs and local equivalences require their separate proofs.",
        },
    }
    output = json.dumps(result, indent=2) + "\n"
    print(output, end="")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
