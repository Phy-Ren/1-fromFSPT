#!/usr/bin/env python3
"""Compute the paper's relative one-form classification in 3+1 dimensions.

Examples:
    python classify.py --orders 4,6 --twist 1,0
    python classify.py --orders 4,6 --twist 1,0 --absolute-bordism
    python classify.py --orders 8 --twist 1
    python classify.py --orders 3,9 --twist 0,0

For A = product_i Z_(N_i), --orders lists N_i and --twist lists the
coefficients of the parity character in the same order. A coefficient is
0 or 1 and must be 0 when its cyclic order is odd. A factor of order 1
is allowed; --orders 1 --twist 0 represents the trivial group.

The default is the paper's relative classification with data (n2, n3, nu4).
The separate --absolute-bordism calculation returns the characters of the
torsion subgroup of absolute twisted-spin bordism. It is a mathematical
comparison, not an established p+ip equivalence of this lattice model.
The old --full spelling is retained as a deprecated alias with a warning.

Successful calls print a JSON object. Invalid arguments produce an error
on standard error and exit with status 2. The Smith presentation and its
input validation are reused from verify_classification_snf.py; this script
and that module must be installed together, with SymPy available.
"""

import argparse
import json
from math import prod
import re
import sys

from verify_classification_snf import presentation, validate_inputs


def integer_csv(value):
    """Parse a nonempty comma-separated list of decimal integers."""
    fields = value.split(",")
    if not fields or any(not field.strip() for field in fields):
        raise argparse.ArgumentTypeError(
            "Expected a nonempty comma-separated list, with no empty entries."
        )
    result = []
    for field in fields:
        token = field.strip()
        if re.fullmatch(r"[+-]?[0-9]+", token) is None:
            raise argparse.ArgumentTypeError(
                f"Expected a decimal integer; received {token!r}."
            )
        try:
            result.append(int(token, 10))
        except ValueError as error:
            raise argparse.ArgumentTypeError(
                "An integer entry exceeds the interpreter's decimal input limit."
            ) from error
    return result


def classify(orders, twist, absolute_bordism=False):
    """Return the relative group or the absolute-bordism comparison.

    The invariant factors describe a product of cyclic groups and omit
    factors of order one. Their coordinates do not specify an embedding
    of the universal Majorana subgroup. The absolute calculation does
    not assert an additional allowed lattice equivalence.
    """
    validate_inputs(orders, twist, absolute_bordism)
    factors = presentation(orders, twist, absolute_bordism)
    nonzero_twist = any(twist)
    if nonzero_twist:
        subgroup_note = (
            "The canonical universal Majorana subgroup has order 16 in the "
            "relative one-form classification. It need not be a direct factor, "
            "and its embedding in these invariant-factor coordinates is not computed. "
            + ("The absolute-bordism comparison takes the mathematical quotient "
               "by this subgroup; no extra lattice equivalence is asserted."
               if absolute_bordism else
               "This subgroup is retained in the paper's classification.")
        )
    else:
        subgroup_note = (
            "For the zero twist the universal Majorana subgroup is trivial. "
            "The relative group and the absolute-bordism torsion-character "
            "comparison have the same invariant factors."
        )
    return {
        "symmetry": {
            "cyclic_orders": list(orders),
            "parity_character": list(twist),
        },
        "invariant_factors": list(factors),
        "group": " x ".join(f"Z_{n}" for n in factors) or "trivial",
        "order": prod(factors),
        "convention": (
            "absolute_bordism_torsion_characters"
            if absolute_bordism else "relative_oneform"
        ),
        "mathematical_scope": {
            "spacetime_dimension": "3+1D",
            "symmetry": "Finite Abelian pure one-form symmetry; all operations are unitary.",
            "twist": "A fixed parity character in Hom(A, Z_2).",
            "classification": (
                "Characters of the torsion subgroup of absolute twisted-spin bordism."
                if absolute_bordism else
                "The paper's relative classification with data (n2, n3, nu4): "
                "Pontryagin dual of twisted bordism modulo the image of ordinary spin bordism."
            ),
            "interpretation": (
                "Mathematical comparison only; an additional symmetric local "
                "one-form p+ip equivalence has not been established."
                if absolute_bordism else
                "The classification problem specified by the paper's one-form lattice construction."
            ),
        },
        "universal_subgroup": {
            "order_in_relative_theory": 16 if nonzero_twist else 1,
            "quotiented_in_this_calculation": bool(absolute_bordism and nonzero_twist),
            "note": subgroup_note,
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Compute the paper's relative 3+1D classification for a finite "
            "Abelian unitary one-form symmetry and a parity character."
        ),
        epilog=(
            "Example: python classify.py --orders 4,6 --twist 1,0. "
            "Use --absolute-bordism for a separate mathematical comparison."
        ),
    )
    parser.add_argument(
        "--orders", type=integer_csv, required=True, metavar="N1,N2,...",
        help="Positive cyclic orders of A, separated by commas.",
    )
    parser.add_argument(
        "--twist", type=integer_csv, required=True, metavar="W1,W2,...",
        help="One coefficient 0 or 1 per factor; coefficients on odd orders must be 0.",
    )
    comparison = parser.add_mutually_exclusive_group()
    comparison.add_argument(
        "--absolute-bordism", action="store_true",
        help="Compute absolute-bordism torsion characters as a mathematical comparison.",
    )
    comparison.add_argument(
        "--full", dest="legacy_full", action="store_true",
        help="Deprecated alias for --absolute-bordism; does not mean a full lattice equivalence.",
    )
    args = parser.parse_args(argv)
    try:
        result = classify(args.orders, args.twist, args.absolute_bordism or args.legacy_full)
    except ValueError as error:
        parser.error(str(error))
    if args.legacy_full:
        print(
            "Warning: --full is deprecated; use --absolute-bordism. "
            "This computes a mathematical bordism comparison, not an established "
            "p+ip equivalence of the one-form lattice model.",
            file=sys.stderr,
        )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
