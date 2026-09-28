#!/usr/bin/env python3
"""Check the 2+1D one-form Kitaev-loop cochain equivalence exactly.

These checks verify the top-cochain correction and its trivialization.
They do not reconstruct the auxiliary-fermion circuit or replace the
cohomological completeness proof. Requires NumPy and the adjacent
cochain_oneform_verify.py module.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from cochain_oneform_verify import cyclic_background, cup, mismatch


def check_loop_equivalence():
    # Each binary 2-cocycle on a 4-simplex has six independent triangle bits.
    ids = np.arange(2**12, dtype=np.int64)
    occupation = cyclic_background(2, 4, ids % 64).reduce()
    twist = cyclic_background(2, 4, ids // 64).reduce()
    shifted = occupation + twist
    obstruction = cup(occupation, occupation) + cup(twist, occupation)
    shifted_obstruction = cup(shifted, shifted) + cup(twist, shifted)
    correction = cup(occupation, twist, 1)
    change = correction.differential() - (shifted_obstruction - obstruction)
    # Applying the equivalence twice leaves an ordinary phase coboundary.
    double_correction = correction + cup(shifted, twist, 1)
    square_boundary = double_correction.lift().scale(2) - twist.lift().differential()

    # The candidate n2=twist is removed, including its top phase, on all tetrahedra.
    tetrahedron_twist = cyclic_background(2, 3, np.arange(8, dtype=np.int64)).reduce()
    trivialization = cup(tetrahedron_twist, tetrahedron_twist, 1).lift().scale(2)
    trivialization = trivialization - tetrahedron_twist.lift().differential()
    return {
        'independent_occupation_twist_pairs': 4096,
        'binary_tetrahedra': 8,
        'failures': {
            'change_of_obstruction': mismatch(change, 2),
            'double_equivalence_is_coboundary': mismatch(square_boundary, 4, top=False),
            'twisted_candidate_top_phase_is_exact': mismatch(trivialization, 4),
        },
        'scope': 'Exact cochain identities; the microscopic equivalence and completeness are proved separately.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check_loop_equivalence()
    output = json.dumps(result, indent=2) + '\n'
    print(output, end='')
    if args.output:
        args.output.write_text(output)
    assert not any(result['failures'].values())
