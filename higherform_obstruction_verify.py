#!/usr/bin/env python3
"""Independent detectors for the higher-form fermion obstructions.

Python standard library only. Steenrod operations are computed in the known
polynomial cohomology of a product of RP-infinity spaces, not by attempting to
infer cohomology from randomly sampled simplices. Separate simplicial checks
verify the integral-Bockstein backgrounds that realize these detectors.

Run: python higherform_obstruction_verify.py --output results/higherform_obstruction_checks.json
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import random


def poly_add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def poly_mul(left, right):
    result = set()
    for a in left:
        for b in right:
            term = tuple(x + y for x, y in zip(a, b))
            result.symmetric_difference_update({term})
    return result


def allocations(total, bounds):
    if not bounds:
        if total == 0:
            yield ()
        return
    for first in range(min(total, bounds[0]) + 1):
        for tail in allocations(total - first, bounds[1:]):
            yield (first,) + tail


def sq(index, polynomial):
    """Sq^index on F2[t1,...,tr], where each ti has degree one."""
    result = set()
    for exponents in polynomial:
        for shifts in allocations(index, exponents):
            if all(math.comb(a, b) % 2 for a, b in zip(exponents, shifts)):
                result.symmetric_difference_update({
                    tuple(a + b for a, b in zip(exponents, shifts))
                })
    return result


def readable(polynomial):
    terms = []
    for exponents in sorted(polynomial, reverse=True):
        factors = [
            f"t{i + 1}" + (f"^{power}" if power != 1 else "")
            for i, power in enumerate(exponents) if power
        ]
        terms.append("*".join(factors) or "1")
    return " + ".join(terms) or "0"


def polynomial_detectors():
    degree_three = sq(1, {(1, 1)})
    degree_four = sq(1, {(1, 1, 1)})
    fundamental_obstruction = sq(1, sq(2, degree_three))
    carry_obstruction = sq(1, sq(2, degree_four))
    expected_four = set(itertools.permutations((4, 2, 1)))
    checks = {
        "degree_three_reduction": degree_three == {(2, 1), (1, 2)},
        "degree_three_Adem": fundamental_obstruction == sq(3, degree_three),
        "degree_three_top_square":
            fundamental_obstruction == poly_mul(degree_three, degree_three),
        "degree_three_nonzero": bool(fundamental_obstruction),
        "degree_four_Adem": carry_obstruction == sq(3, degree_four),
        "degree_four_six_distinct_monomials": carry_obstruction == expected_four,
        "degree_four_nonzero": bool(carry_obstruction),
    }
    return {
        "rho2_beta2_t1t2": readable(degree_three),
        "Sq1_Sq2_rho2_beta2_t1t2": readable(fundamental_obstruction),
        "rho2_beta2_t1t2t3": readable(degree_four),
        "Sq1_Sq2_rho2_beta2_t1t2t3": readable(carry_obstruction),
        "checks": checks,
    }


def faces(dim, degree):
    return tuple(itertools.combinations(range(dim + 1), degree + 1))


def coboundary(values, dim, degree):
    return {
        face: sum(
            (-1)**j * values[face[:j] + face[j + 1:]]
            for j in range(len(face))
        )
        for face in faces(dim, degree + 1)
    }


def product_of_degree_one_classes(vertices, degree):
    """Alexander-Whitney product t1 cup ... cup t_degree."""
    result = {}
    for face in faces(len(vertices) - 1, degree):
        value = 1
        for bit in range(degree):
            value *= (
                (vertices[face[bit]] ^ vertices[face[bit + 1]]) >> bit
            ) & 1
        result[face] = value
    return result


def integral_bockstein(values, dim, degree):
    derivative = coboundary(values, dim, degree)
    assert all(value % 2 == 0 for value in derivative.values())
    return {face: value // 2 for face, value in derivative.items()}


def simplicial_checks():
    seed = 20260928
    rng = random.Random(seed)
    samples = 256
    orders = (2, 4, 6, 8, 12, 16)
    failures = {
        "integral_degree_three_is_closed": 0,
        "integral_degree_four_is_closed": 0,
        "fundamental_degree_three_background_is_closed": 0,
        "fundamental_degree_three_parity_is_correct": 0,
        "carry_background_is_closed": 0,
        "carry_is_exactly_beta2_t1t2t3": 0,
        "fundamental_degree_four_background_is_closed": 0,
        "fundamental_degree_four_parity_is_correct": 0,
    }
    dim = 7
    for _ in range(samples):
        vertices = (0,) + tuple(rng.randrange(8) for _ in range(dim))
        p2 = product_of_degree_one_classes(vertices, 2)
        p3 = product_of_degree_one_classes(vertices, 3)
        z3 = integral_bockstein(p2, dim, 2)
        z4 = integral_bockstein(p3, dim, 3)
        failures["integral_degree_three_is_closed"] += sum(
            bool(v) for v in coboundary(z3, dim, 3).values()
        )
        failures["integral_degree_four_is_closed"] += sum(
            bool(v) for v in coboundary(z4, dim, 4).values()
        )
        for order in orders:
            a3 = {f: v % order for f, v in z3.items()}
            a4 = {f: v % order for f, v in z4.items()}
            carry_a3 = {f: (order // 2) * v for f, v in p3.items()}
            carry_derivative = coboundary(carry_a3, dim, 3)
            failures["fundamental_degree_three_background_is_closed"] += sum(
                bool(v % order) for v in coboundary(a3, dim, 3).values()
            )
            failures["fundamental_degree_four_background_is_closed"] += sum(
                bool(v % order) for v in coboundary(a4, dim, 4).values()
            )
            failures["carry_background_is_closed"] += sum(
                bool(v % order) for v in carry_derivative.values()
            )
            failures["fundamental_degree_three_parity_is_correct"] += sum(
                (a3[f] - v) % 2 != 0 for f, v in z3.items()
            )
            failures["fundamental_degree_four_parity_is_correct"] += sum(
                (a4[f] - v) % 2 != 0 for f, v in z4.items()
            )
            failures["carry_is_exactly_beta2_t1t2t3"] += sum(
                carry_derivative[f] != order * v for f, v in z4.items()
            )
    return {
        "seed": seed, "seven_simplices": samples,
        "cyclic_orders": orders, "failures": failures,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    algebra = polynomial_detectors()
    cochains = simplicial_checks()
    failed = (
        sum(not value for value in algebra["checks"].values())
        + sum(cochains["failures"].values())
    )
    result = {
        "scope": "cohomology detectors and realizing Bockstein backgrounds",
        "polynomial_detectors": algebra,
        "simplicial_background_checks": cochains,
        "total_failures": failed,
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    raise SystemExit(bool(failed))


if __name__ == "__main__":
    main()
