#!/usr/bin/env python3
"""Classify finite-Abelian unitary pure one-form responses in 4+1D.

Examples:
    python classify_41d.py --orders 2 --twist 1
    python classify_41d.py --orders 4,8 --twist 1,0
    python classify_41d.py --orders 6,10 --twist 1,1

The convention is characters of twisted spin bordism in degree five,
trivial on the trivial-background ordinary-spin image. Every class has a
bosonic representative. This calculates an abstract finite group, not a
microscopic circuit or a natural choice of invariant-factor generators.

Only SymPy is required. The presentation is the integral H^6 quotient by
beta_2(Sq^2 + w cup) on H^3(-,Z_2). No 3+1D modules are imported.
"""

import argparse
from itertools import combinations
import json
from math import gcd, prod
import re

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def validate_inputs(orders, twist):
    """Validate a cyclic decomposition and its parity character."""
    if len(orders) != len(twist):
        raise ValueError("Cyclic orders and twist coefficients must have equal lengths.")
    if any(type(n) is not int or n < 1 for n in orders):
        raise ValueError("Every cyclic order must be a positive integer.")
    if any(type(w) is not int or w not in (0, 1) for w in twist):
        raise ValueError("Every twist coefficient must be the integer 0 or 1.")
    if any(n % 2 and w for n, w in zip(orders, twist)):
        raise ValueError("A parity character must vanish on every odd-order factor.")


def presentation(orders, twist):
    """Return Smith factors of the full bosonic quotient presentation.

    D_i = gamma_i^2 has order two for even N_i.
    C_ij = gamma_i gamma_j, i < j, has order gcd(N_i,N_j).
    Each even i contributes the relation
      gamma_i^2 + sum_j w_j (N_j/2) gamma_j gamma_i = 0.
    Integral degree-three generators anticommute; their squares need not
    vanish, and have order two. The oriented mixed coefficients below
    preserve this sign convention before taking the quotient.
    """
    validate_inputs(orders, twist)
    even = [i for i, n in enumerate(orders) if n % 2 == 0]
    diagonal = {i: k for k, i in enumerate(even)}
    pairs = list(combinations(range(len(orders)), 2))
    mixed = {pair: len(even) + k for k, pair in enumerate(pairs)}
    generator_orders = [2] * len(even)
    generator_orders += [gcd(orders[i], orders[j]) for i, j in pairs]
    size = len(generator_orders)
    if not size:
        return ()
    relations = []
    for k, n in enumerate(generator_orders):
        row = [0] * size
        row[k] = n
        relations.append(row)
    for i in even:
        row = [0] * size
        row[diagonal[i]] = 1 + twist[i] * (orders[i] // 2)
        for j, w in enumerate(twist):
            if not w or j == i:
                continue
            pair = tuple(sorted((i, j)))
            sign = 1 if j < i else -1
            row[mixed[pair]] += sign * (orders[j] // 2)
        relations.append(row)
    reduced = smith_normal_form(Matrix(relations), domain=ZZ)
    factors = tuple(abs(int(reduced[k, k])) for k in range(size))
    if any(n == 0 for n in factors):
        raise ArithmeticError("The finite presentation unexpectedly has a free factor.")
    return tuple(n for n in factors if n != 1)


def classify(orders, twist):
    """Describe the abstract response group and the scope of the calculation."""
    factors = presentation(orders, twist)
    extra = any(w and n % 4 == 2 for n, w in zip(orders, twist))
    return {
        "spacetime_dimension": "4+1D",
        "symmetry": {"cyclic_orders": list(orders), "parity_character": list(twist)},
        "invariant_factors": list(factors),
        "group": " x ".join(f"Z_{n}" for n in factors) or "trivial",
        "order": prod(factors),
        "convention": "closed_response_relative_to_trivial_background",
        "intrinsically_fermionic_classes": False,
        "extra_bosonic_order_two_class": extra,
        "mathematical_scope": {
            "symmetry": "Finite Abelian, unitary, pure one-form symmetry.",
            "classification": "Characters of degree-five twisted spin bordism, trivial on the ordinary-spin image.",
            "representatives": "Every class has n2=n3=n4=0 and an ordinary bosonic degree-five response.",
            "integer_layer": "The integral p+ip occupation H^2(X,Z) and its degree-one residual redefinitions H^1(X,Z) vanish; the full twisted-spin AHSS also excludes extra degree-five contributions.",
            "presentation": "H^6(B^2 A,Z) modulo beta_2(Sq^2 kappa + w cup kappa) for kappa in H^3(B^2 A,Z_2).",
            "extra_class": "The bosonic response (1/2) w cup Sq^1 w survives precisely when Sq^1 w is nonzero.",
            "coordinates": "Invariant-factor coordinates are abstract; no natural generator embedding is claimed.",
            "limits": "No continuous symmetry, nontrivial two-group data, or microscopic circuit verification is included.",
        },
    }


def integer_csv(value):
    """Parse a nonempty comma-separated list of decimal integers."""
    fields = value.split(",")
    if any(not field.strip() for field in fields):
        raise argparse.ArgumentTypeError("Expected a nonempty list with no empty entries.")
    output = []
    for field in fields:
        token = field.strip()
        if re.fullmatch(r"[+-]?[0-9]+", token) is None:
            raise argparse.ArgumentTypeError(f"Expected a decimal integer; received {token!r}.")
        try:
            output.append(int(token, 10))
        except ValueError as error:
            raise argparse.ArgumentTypeError("An integer exceeds the interpreter's decimal input limit.") from error
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=integer_csv, required=True, metavar="N1,N2,...")
    parser.add_argument("--twist", type=integer_csv, required=True, metavar="W1,W2,...")
    args = parser.parse_args(argv)
    try:
        result = classify(args.orders, args.twist)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
