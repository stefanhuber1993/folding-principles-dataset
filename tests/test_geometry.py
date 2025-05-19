# tests/test_geometry.py
import sys, pathlib, numpy as np
from math import cos, sin, radians, degrees, acos
from Bio.PDB import Atom, Residue, Vector

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from dataset.geometry import get_cb_from_residue, CA_CB_BOND


def make_ala_residue():
    """Realistic Ala backbone: N-CA 1.458 Å, C-CA 1.525 Å, ∠N-CA-C 110°."""
    res = Residue.Residue((" ", 1, " "), "ALA", "")
    CA = np.array([0.0, 0.0, 0.0])
    N  = np.array([1.458, 0.0, 0.0])
    C  = np.array([
        1.525 * cos(radians(110.0)),
        1.525 * sin(radians(110.0)),
        0.0
    ])
    res.add(Atom.Atom("N",  N,  0, 0, " ", "N",  1))
    res.add(Atom.Atom("CA", CA, 0, 0, " ", "CA", 2))
    res.add(Atom.Atom("C",  C,  0, 0, " ", "C",  3))
    return res


def angle(v1: Vector, v2: Vector) -> float:
    return degrees(acos((v1 * v2) / (v1.norm() * v2.norm())))


def test_pseudo_cb_geometry():
    res = make_ala_residue()
    cb = get_cb_from_residue(res)
    assert cb is not None, "function returned None"

    CA = Vector(res["CA"].get_coord())
    N  = Vector(res["N"].get_coord())
    C  = Vector(res["C"].get_coord())

    # 1 · bond length
    assert np.isclose((cb - CA).norm(), CA_CB_BOND, atol=1e-3)

    # 2 · the two angles must be equal and in tetrahedral range
    theta_N = angle(N - CA, cb - CA)
    theta_C = angle(C - CA, cb - CA)
    assert abs(theta_N - theta_C) < 0.5           # equality
    assert 95.0 < theta_N < 115.0                 # reasonable tetrahedral

    # 3 · Cβ on side-chain half-space (opposite carbonyl C)
    assert (cb - CA) * (C - CA) < 0
