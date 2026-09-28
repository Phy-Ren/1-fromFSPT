#!/usr/bin/env python3
"""Exact 4+1D finite one-form coboundary checks.

Requires Python 3.10+, NumPy, and the companion cochain_oneform_verify.py.
Place both files in the same directory, or put the companion directory on
PYTHONPATH. All computations use integer numerators modulo 4 or 16.

Run: python cochain_41d_verify.py --output cochain_41d_checks.json

This script checks cochain identities; the completeness proof and the finite
Abelian classification are separate. It does not regenerate lattice projector
data or claim exhaustive enumeration of all finite groups.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from cochain_oneform_verify import (
    cup, cyclic_background, mismatch, phase, pontryagin, word,
)


def adem_x2(x):
    terms = [word(w, x, x, x, x) for w in
             ("1213243", "1213431", "1232141", "1234321")]
    return terms[0] + terms[1] + terms[2] + terms[3]


def majorana_boundary(x, omega):
    """Return q=x^2+omega*x and four times O5^gamma[x;omega]."""
    x2 = cup(x, x)
    omega_x = cup(omega, x)
    q = x2 + omega_x
    B = x.beta()
    half_part = (
        word("1231343", omega, omega, x, x)
        + adem_x2(x)
        + cup(x2, omega_x, 3)
    )
    quarter_part = (
        cup(omega.lift(), B, integer=True)
        + cup(B, B, 1, integer=True)
    )
    return q, phase((2, half_part), (1, quarter_part), mod=4)


def obstruction(q, omega):
    """Return four times the additive complex-fermion obstruction."""
    return phase((2, cup(q, q, 2) + cup(omega, q)), mod=4)


def check_boundary(x, omega):
    q, gamma = majorana_boundary(x, omega)
    return {
        "closed_complex_occupation": mismatch(q.differential(), 2),
        "full_majorana_coboundary": mismatch(
            gamma.differential() - obstruction(q, omega), 4
        ),
    }


def run_checks():
    cases = {}
    ids = np.arange(2**15, dtype=np.int64)
    x = cyclic_background(2, 6, ids).reduce()
    zero = x.scale(0).reduce()
    cases["all_binary_six_simplices_split"] = {
        "simplices": len(ids), "failures": check_boundary(x, zero),
    }
    cases["all_binary_six_simplices_diagonal_twist"] = {
        "simplices": len(ids), "failures": check_boundary(x, x),
    }

    seed = 20260928
    rng = np.random.default_rng(seed)
    count = 8192
    x = cyclic_background(2, 6, rng.integers(2**15, size=count)).reduce()
    omega = cyclic_background(2, 6, rng.integers(2**15, size=count)).reduce()
    cases["independent_binary_source_and_twist"] = {
        "simplices": count, "seed": seed,
        "failures": check_boundary(x, omega),
    }

    y = cyclic_background(2, 6, rng.integers(2**15, size=count)).reduce()
    q, gamma = majorana_boundary(x, omega)
    q_other, gamma_other = majorana_boundary(y, omega)
    stacked_phase = phase(
        (1, gamma), (1, gamma_other), (2, cup(q, q_other, 3)), mod=4
    )
    cases["complex_stacking_phase"] = {
        "simplices": count,
        "failures": {
            "stacking_preserves_phase_equation": mismatch(
                stacked_phase.differential()
                - obstruction(q + q_other, omega), 4
            ),
        },
    }

    # The sole possible secondary source is x=omega. Its phase is already
    # integrated by the universal 3+1D Majorana root in the same convention.
    ids = np.arange(2**10, dtype=np.int64)
    u = cyclic_background(2, 5, ids).reduce()
    q, gamma = majorana_boundary(u, u)
    b = cup(u, u, 1)
    root_phase = phase(
        (3, pontryagin(u)), (8, cup(b, u, 1)), mod=16
    )
    cases["universal_secondary_source"] = {
        "simplices": len(ids),
        "failures": {
            "vanishing_complex_occupation": mismatch(q, 2, top=False),
            "secondary_phase_is_exact": mismatch(
                root_phase.differential() - gamma.scale(4), 16
            ),
        },
    }
    return {
        "convention": "four-word x2 and signed integer interval-cut cups",
        "dependency": "companion cochain_oneform_verify.py",
        "scope": "cochain identities; no regenerated lattice projector data",
        "cases": cases,
        "total_failures": sum(
            sum(case["failures"].values()) for case in cases.values()
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")
    raise SystemExit(bool(result["total_failures"]))


if __name__ == "__main__":
    main()
