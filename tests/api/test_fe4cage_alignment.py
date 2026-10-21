"""Regression test for the equal-weight axis frame in ``min_rmsd``.

``tests/data/Fe4Cage.xyz`` is a T-symmetric Fe4L6 cage (316 atoms). No atom
has a unique canonical rank, so ``min_rmsd`` aligns the structures on their
principal axes and searches the axis flips. T symmetry makes all three
principal moments equal (with mass weights and with equal weights alike), so
the principal axes are an arbitrary choice and need not correspond between the
two structures. With mass-weighted axes alone, one of the 11 non-identity T
operations (a C3) scored 2.007 A instead of 0.0 A. The equal-weight frame gives
the search a second, independent set of axes and recovers it.

The point-group operations are generated at test time with
``get_symmetry_operations``, an alignment-independent geometric symmetry
search. Each one is a permutation that maps the cage exactly onto itself, so
``get_irmsd`` must return ~0 A for every one of them.
"""

from pathlib import Path

import pytest

from irmsd.api.irmsd_exposed import get_irmsd
from irmsd.api.symmetry_exposed import get_symmetry_operations
from irmsd.utils.xyz import read_extxyz

FE4CAGE_XYZ = Path(__file__).parents[1] / "data" / "Fe4Cage.xyz"


@pytest.fixture(scope="module")
def fe4cage():
    mol = read_extxyz(FE4CAGE_XYZ)
    return mol.get_atomic_numbers(), mol.get_positions()


@pytest.fixture(scope="module")
def fe4cage_operations(fe4cage):
    Z, pos = fe4cage
    symbol, ops = get_symmetry_operations(Z, pos, max_atoms=None)
    return symbol, [op for op in ops if op.label != "E"]


def test_fe4cage_point_group(fe4cage_operations):
    symbol, ops = fe4cage_operations
    assert symbol == "T"
    assert sorted(op.label for op in ops) == sorted(
        ["C2"] * 3 + ["C3"] * 4 + ["C3^2"] * 4
    )


def test_fe4cage_pointgroup_isomers_score_zero(fe4cage, fe4cage_operations):
    Z, pos = fe4cage
    _, ops = fe4cage_operations
    rmsds = {}
    for i, op in enumerate(ops):
        rmsdval, *_ = get_irmsd(Z, pos, Z, pos[op.permutation])
        rmsds[f"{i}:{op.label}"] = float(rmsdval)

    failed = {k: v for k, v in rmsds.items() if v >= 1e-4}
    assert not failed, (
        f"{len(failed)}/{len(ops)} exact point-group operations did not score "
        f"~0 A: {failed} -- the equal-weight axis frame may be missing or broken"
    )
