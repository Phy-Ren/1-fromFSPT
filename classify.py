#!/usr/bin/env python3
"""Classify finite Abelian unitary one-form symmetries in 3+1 dimensions.

Examples:
    python classify.py --orders 4,6 --twist 1,0
    python classify.py --orders 4,6 --twist 1,0 --full
    python classify.py --orders 8 --twist 1
    python classify.py --orders 3,9 --twist 0,0

For A = product_i Z_(N_i), --orders lists N_i and --twist lists the
coefficients of the parity character in the same order. A coefficient is
0 or 1 and must be 0 when its cyclic order is odd. A factor of order 1
is allowed; --orders 1 --twist 0 represents the trivial group.

The default is the reduced three-layer classification. The --full option
also includes the integer p+ip identifications, equivalently quotienting
by the canonical universal order-sixteen subgroup for a nonzero twist.
These are different equivalence conventions, including for Z_2 and Z_4.

Successful calls print a JSON object. Invalid arguments produce an error
on standard error and exit with status 2. The Smith presentation and its
input validation are reused from verify_classification_snf.py; this script
and that module must be installed together, with SymPy available.
"""

import argparse
import json
from math import prod
import re

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


def classify(orders, twist, full=False):
    """Return the abstract classification and its equivalence convention.

    The invariant factors describe a product of cyclic groups and omit
    factors of order one. Their coordinates do not specify an embedding
    of the universal Majorana subgroup.
    """
    validate_inputs(orders, twist, full)
    factors = presentation(orders, twist, full)
    nonzero_twist = any(twist)
    if nonzero_twist:
        subgroup_note = (
            "The canonical pullback of the universal Z_2 Majorana group has "
            "order 16 in the reduced theory. It need not be a direct factor. "
            "Its embedding in these invariant-factor coordinates is not computed. "
            + ("This subgroup has been quotiented out." if full else
               "This subgroup is retained; --full quotients it out.")
        )
    else:
        subgroup_note = (
            "For the zero twist the universal Majorana subgroup is trivial. "
            "Including integer-layer identifications leaves this group unchanged."
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
            "full_integer_layer_quotient" if full else "three_layer_reduced"
        ),
        "mathematical_scope": {
            "spacetime_dimension": "3+1D",
            "symmetry": "Finite Abelian pure one-form symmetry; all operations are unitary.",
            "twist": "A fixed parity character in Hom(A, Z_2).",
            "classification": (
                "The torsion deformation group after integer p+ip identifications."
                if full else
                "The three-layer group with integer p+ip redefinitions held fixed; "
                "Pontryagin dual of twisted bordism modulo the image of ordinary spin bordism."
            ),
        },
        "universal_subgroup": {
            "order_in_reduced_theory": 16 if nonzero_twist else 1,
            "quotiented_out": bool(full and nonzero_twist),
            "note": subgroup_note,
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Compute the 3+1D classification for a finite Abelian unitary "
            "one-form symmetry and a specified parity character."
        ),
        epilog=(
            "Example: python classify.py --orders 4,6 --twist 1,0. "
            "Append --full to include the integer p+ip identifications."
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
    parser.add_argument(
        "--full", action="store_true",
        help="Also quotient by the integer p+ip identifications.",
    )
    args = parser.parse_args(argv)
    try:
        result = classify(args.orders, args.twist, args.full)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
