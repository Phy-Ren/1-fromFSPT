"""Check the one-form classification by two independent group presentations.

The first construction uses all quadratic generators in the original cyclic
coordinates, the Wu relations, and the universal Majorana extension. The
second separates primary factors and uses the adapted-character product
formula. Both the paper's relative one-form classification and the mathematical
absolute-bordism torsion-character comparison are checked. The latter is not
an established p+ip equivalence of the lattice model. These arithmetic checks
do not establish the topological inputs or microscopic equivalences.
"""

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import platform

import sympy
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def validate_inputs(ns, ws, absolute_bordism=False):
    """Validate cyclic orders and a well-defined parity character."""
    if len(ns) != len(ws):
        raise ValueError("Cyclic orders and character coefficients must have equal lengths.")
    if any(type(n) is not int or n < 1 for n in ns):
        raise ValueError("Every cyclic order must be a positive integer.")
    if any(type(w) is not int or w not in (0, 1) for w in ws):
        raise ValueError("Every character coefficient must be the integer 0 or 1.")
    if any(n % 2 and w for n, w in zip(ns, ws)):
        raise ValueError("A parity character must vanish on every odd-order factor.")
    if type(absolute_bordism) is not bool:
        raise ValueError("The absolute-bordism comparison flag must be a boolean.")


def presentation(ns, ws, absolute_bordism=False):
    """Return Smith factors for the relative group or its absolute comparison.

    ``absolute_bordism=True`` computes characters of absolute bordism torsion.
    It does not assert an extra physical equivalence of the one-form model.
    """
    validate_inputs(ns, ws, absolute_bordism)
    r = len(ns)
    pairs = list(combinations(range(r), 2))
    mixed = {pair: r + k for k, pair in enumerate(pairs)}
    orders = [2 * n if n % 2 == 0 else n for n in ns]
    orders += [gcd(ns[i], ns[j]) for i, j in pairs]
    nonzero = any(ws)
    size = len(orders) + int(nonzero)
    if size == 0:
        return ()
    rows = []
    for i, n in enumerate(orders):
        row = [0] * size
        row[i] = n
        rows.append(row)
    for i, n in enumerate(ns):
        if n % 2:
            continue
        row = [0] * size
        # 1 - w_i and 1 + w_i agree modulo the order 2N_i.
        row[i] = n * (1 - ws[i])
        for j in range(r):
            if i != j and ws[j]:
                row[mixed[tuple(sorted((i, j)))]] = gcd(n, ns[j]) // 2
        rows.append(row)
    if nonzero:
        complex_layer = any(w and n % 4 == 2 for n, w in zip(ns, ws))
        denominator = 2 if complex_layer else 4
        row = [0] * size
        row[-1] = 4 if complex_layer else 2
        for i, n in enumerate(ns):
            if ws[i]:
                row[i] = -(n // denominator)
        for (i, j), k in mixed.items():
            if ws[i] and ws[j]:
                row[k] = -(gcd(ns[i], ns[j]) // denominator)
        rows.append(row)
        if absolute_bordism:
            row = [0] * size
            row[-1] = 1
            rows.append(row)
    snf = smith_normal_form(Matrix(rows), domain=ZZ)
    return tuple(
        abs(int(snf[i, i]))
        for i in range(size)
        if abs(snf[i, i]) != 1
    )


def two_primary_formula(exponents, ws, absolute_bordism=False):
    """Return invariant factors for two-primary cyclic factors."""
    validate_inputs([2 ** e for e in exponents], ws, absolute_bordism)
    if any(type(e) is not int or e < 1 for e in exponents):
        raise ValueError("Two-primary exponents must be positive integers.")
    if not any(ws):
        factors = list(exponents)
        factors += [min(a, b) for a, b in combinations(exponents, 2)]
        return tuple(2 ** e for e in sorted(factors))
    i0 = min((i for i, w in enumerate(ws) if w), key=lambda i: exponents[i])
    k = exponents[i0]
    spectators = [e for i, e in enumerate(exponents) if i != i0]
    if absolute_bordism:
        factors = [] if k <= 2 else [k - 2]
    else:
        factors = [4] if k <= 2 else [1, k + 1]
    for s in spectators:
        factors += [s + 1, min(k, s) - 1]
    factors += [min(a, b) for a, b in combinations(spectators, 2)]
    return tuple(2 ** e for e in sorted(factors) if e)


def invariant_factors(orders):
    """Combine primary powers without calling Smith normal form."""
    prime_powers = {}
    for n in orders:
        for prime, exponent in sympy.factorint(n).items():
            prime_powers.setdefault(int(prime), []).append(int(exponent))
    length = max((len(exponents) for exponents in prime_powers.values()), default=0)
    result = [1] * length
    for prime, exponents in prime_powers.items():
        aligned = [0] * (length - len(exponents)) + sorted(exponents)
        for i, exponent in enumerate(aligned):
            result[i] *= prime ** exponent
    return tuple(result)


def primary_formula(ns, ws, absolute_bordism=False):
    """Use coprime splitting, odd quadratic factors, and an adapted twist."""
    validate_inputs(ns, ws, absolute_bordism)
    odd_orders = []
    exponents = []
    two_character = []
    for n, w in zip(ns, ws):
        exponent = 0
        odd = n
        while odd % 2 == 0:
            odd //= 2
            exponent += 1
        odd_orders.append(odd)
        if exponent:
            exponents.append(exponent)
            two_character.append(w)
    factors = odd_orders + [gcd(a, b) for a, b in combinations(odd_orders, 2)]
    factors += list(two_primary_formula(exponents, two_character, absolute_bordism))
    return invariant_factors(factors)


def check_case(ns, ws, absolute_bordism):
    actual = presentation(ns, ws, absolute_bordism)
    expected = primary_formula(ns, ws, absolute_bordism)
    if actual != expected:
        raise AssertionError({
            "cyclic_orders": ns,
            "parity_character": ws,
            "absolute_bordism_comparison": absolute_bordism,
            "presentation": actual,
            "primary_formula": expected,
        })


def valid_characters(ns):
    choices = [(0,) if n % 2 else (0, 1) for n in ns]
    return product(*choices)


def run_checks():
    counts = {"two_primary_presentations": 0, "cyclic_presentations": 0,
              "odd_and_mixed_presentations": 0, "trivial_presentations": 0,
              "invalid_inputs_rejected": 0}
    for length in range(1, 4):
        for exponents in product(range(1, 5), repeat=length):
            ns = [2 ** e for e in exponents]
            for ws in valid_characters(ns):
                for absolute_bordism in (False, True):
                    check_case(ns, ws, absolute_bordism)
                    counts["two_primary_presentations"] += 1
    for n in range(1, 33):
        for ws in valid_characters([n]):
            for absolute_bordism in (False, True):
                check_case([n], ws, absolute_bordism)
                counts["cyclic_presentations"] += 1
        if n % 2 == 0:
            reduced = ((8 * n,) if n % 4 == 2 else
                       (4 * n,) if n % 8 == 4 else (2, 2 * n))
            order = n // 2 if n % 4 == 2 else n // 4
            absolute_group = () if order == 1 else (order,)
            assert presentation([n], [1]) == reduced
            assert presentation([n], [1], True) == absolute_group
    mixed_orders = (2, 3, 4, 6, 9, 10, 12, 15)
    for length in range(1, 4):
        for ns in product(mixed_orders, repeat=length):
            # Exercise an odd-primary sector in every case of this family.
            if all(n & (n - 1) == 0 for n in ns):
                continue
            for ws in valid_characters(ns):
                for absolute_bordism in (False, True):
                    check_case(ns, ws, absolute_bordism)
                    counts["odd_and_mixed_presentations"] += 1
    for absolute_bordism in (False, True):
        check_case([], [], absolute_bordism)
        check_case([1, 1], [0, 0], absolute_bordism)
        counts["trivial_presentations"] += 2
    invalid = [([3], [1], False), ([2, 3], [1, 1], False),
               ([2], [], False), ([0], [0], False), ([-2], [0], False),
               ([2.0], [0], False), ([2], [2], False), ([2], [0], 1)]
    for ns, ws, absolute_bordism in invalid:
        try:
            presentation(ns, ws, absolute_bordism)
        except ValueError:
            counts["invalid_inputs_rejected"] += 1
        else:
            raise AssertionError(f"Invalid input was accepted: {ns}, {ws}, {absolute_bordism}")
    return {
        "status": "passed",
        "failures": 0,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "sympy_version": sympy.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "counts": counts,
        "coverage": {
            "two_primary": "All characters on one to three cyclic factors with exponents 1 through 4.",
            "cyclic": "All valid characters at orders 1 through 32, with an independent cyclic closed-form check.",
            "odd_and_mixed": "All valid characters on one to three factors of orders 2, 3, 4, 6, 9, 10, 12, 15 that contain an odd-primary part.",
            "calculations": ["relative one-form classification", "absolute-bordism torsion-character comparison"],
            "interpretation": "The absolute comparison is mathematical; no additional one-form lattice p+ip equivalence is asserted.",
            "limitation": "Checks integer presentations and invariant factors; does not replace the topological completeness proof.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the provenance summary as JSON.")
    args = parser.parse_args()
    summary = run_checks()
    counts = summary["counts"]
    print(f"Two-primary presentations checked: {counts['two_primary_presentations']}.")
    print(f"Cyclic presentations checked: {counts['cyclic_presentations']}.")
    print(f"Odd and mixed-order presentations checked: {counts['odd_and_mixed_presentations']}.")
    print(f"Trivial-group presentations checked: {counts['trivial_presentations']}.")
    print(f"Invalid inputs rejected: {counts['invalid_inputs_rejected']}.")
    print("All checks passed; zero failures.")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"Provenance summary written to {args.output}.")


if __name__ == "__main__":
    main()
