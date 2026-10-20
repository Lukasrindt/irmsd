"""Regression test for the equal-mass axis fallback in ``min_rmsd``.

``tests/data/Fe4Cage.xyz`` is a T-symmetric Fe4L6 cage (316 atoms): 4 heavy Fe
centers among light organic linkers (C/H/N/O). No atom has a unique canonical
rank, so ``min_rmsd`` falls back to principal-axis pre-alignment plus the
32-way grid search. With mass-weighted axes, the 4 Fe atoms dominate the
inertia tensor and skew the axes, so before the fallback one of the 11
non-identity T operations (a C3) scored 2.007 A instead of 0.0 A. This is
unrelated to the degenerate-axis limitation of high-symmetry molecules
(C60, adamantane): the cage's axes are well defined, only biased by mass.

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
        f"~0 A: {failed} -- the equal-mass axis fallback may be missing or broken"
    )
