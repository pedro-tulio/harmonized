"""Subspace harmonisation utilities for adaptive MDeNM/pyAdMD.

The functions operate only on the atom subset used to define the normal-mode
subspace (typically protein C-alpha atoms). They do not alter the MD system.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class SubspaceResult:
    q_new: np.ndarray
    harmonized_modes: np.ndarray
    singular_values: np.ndarray
    principal_angles_deg: np.ndarray
    omega: float
    rotation: np.ndarray
    matching_coefficients: np.ndarray


def _orthonormal_columns(modes: np.ndarray) -> np.ndarray:
    """Return an orthonormal column basis from modes shaped (k, n_atoms, 3)."""
    modes = np.asarray(modes, dtype=float)
    if modes.ndim != 3 or modes.shape[2] != 3:
        raise ValueError("modes must have shape (k, n_atoms, 3)")
    a = modes.reshape(modes.shape[0], -1).T
    q, _ = np.linalg.qr(a, mode="reduced")
    return q


def kabsch_rotation(mobile: np.ndarray,
                    reference: np.ndarray,
                    weights: np.ndarray | None = None) -> np.ndarray:
    """Rotation R such that centered mobile @ R best matches centered reference."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)

    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("mobile/reference must both have shape (n_atoms, 3)")

    if weights is None:
        weights = np.ones(mobile.shape[0], dtype=float)

    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or len(weights) != mobile.shape[0]:
        raise ValueError("weights must contain one value per selected atom")
    if np.any(weights < 0.0) or weights.sum() <= 0.0:
        raise ValueError("weights must be non-negative and have a positive sum")

    weights = weights / weights.sum()

    cm = np.sum(mobile * weights[:, None], axis=0)
    cr = np.sum(reference * weights[:, None], axis=0)

    x = mobile - cm
    y = reference - cr

    cov = (x * weights[:, None]).T @ y
    u, _, vt = np.linalg.svd(cov)

    r = u @ vt
    if np.linalg.det(r) < 0.0:
        u[:, -1] *= -1.0
        r = u @ vt

    return r


def harmonize_and_transport(old_modes: np.ndarray,
                            new_modes: np.ndarray,
                            q_old: np.ndarray,
                            old_coords: np.ndarray,
                            new_coords: np.ndarray,
                            masses: np.ndarray | None = None,
                            alpha: float = 1.0) -> SubspaceResult:
    """Track an active old subspace U_m inside a possibly larger new subspace V_M.

    ``old_modes`` defines the active subspace of dimension m.
    ``new_modes`` defines the new ENM search subspace of dimension M, with M >= m.

    The rectangular overlap matrix U_m^T V_M is decomposed by thin SVD. Its
    right singular vectors select the m-dimensional part of V_M with maximal
    overlap with U_m. That matched subspace is then Procrustes-harmonized and
    used to transport the current physical excitation direction Q.
    """
    old_modes = np.asarray(old_modes, dtype=float)
    new_modes = np.asarray(new_modes, dtype=float)
    q_old = np.asarray(q_old, dtype=float)
    old_coords = np.asarray(old_coords, dtype=float)
    new_coords = np.asarray(new_coords, dtype=float)

    if old_modes.ndim != 3 or old_modes.shape[-1] != 3:
        raise ValueError("old_modes must have shape (m, n_atoms, 3)")
    if new_modes.ndim != 3 or new_modes.shape[-1] != 3:
        raise ValueError("new_modes must have shape (M, n_atoms, 3)")
    if old_modes.shape[1:] != new_modes.shape[1:]:
        raise ValueError(
            "old_modes and new_modes must contain the same atoms in the same order"
        )

    m = old_modes.shape[0]
    M = new_modes.shape[0]

    if M < m:
        raise ValueError(f"search subspace must satisfy M >= m; got M={M}, m={m}")
    if q_old.shape != old_modes.shape[1:]:
        raise ValueError("q_old must have shape (n_atoms, 3)")
    if old_coords.shape != q_old.shape or new_coords.shape != q_old.shape:
        raise ValueError("old_coords/new_coords must have shape (n_atoms, 3)")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be between 0 and 1")

    r = kabsch_rotation(new_coords, old_coords, masses)
    new_modes_rot = new_modes @ r

    u = _orthonormal_columns(old_modes)   # (3N, m)
    v = _orthonormal_columns(new_modes_rot)  # (3N, M)

    s = u.T @ v                           # (m, M)

    l, sigma, rt = np.linalg.svd(s, full_matrices=False)
    r_m = rt.T                            # (M, m)

    v_match = v @ r_m                     # (3N, m)
    v_harm = v_match @ l.T                # (3N, m)

    q0 = q_old.reshape(-1)
    q0_norm = np.linalg.norm(q0)
    if q0_norm < 1e-14:
        raise ValueError("q_old has zero norm")
    q0 /= q0_norm

    q_transport = v_harm @ (v_harm.T @ q0)
    nt = np.linalg.norm(q_transport)
    if nt < 1e-12:
        raise RuntimeError(
            "transported Q is nearly zero: search subspace has insufficient overlap"
        )
    q_transport /= nt

    if np.dot(q_transport, q0) < 0.0:
        q_transport *= -1.0

    q_mix = (1.0 - alpha) * q0 + alpha * q_transport
    nq = np.linalg.norm(q_mix)
    if nq < 1e-12:
        q_mix = q_transport
        nq = 1.0

    q_new = (q_mix / nq).reshape(q_old.shape)

    sigma_clip = np.clip(sigma, 0.0, 1.0)
    angles = np.degrees(np.arccos(sigma_clip))
    omega = float(np.mean(sigma_clip ** 2))

    modes_harm = v_harm.T.reshape((m,) + old_modes.shape[1:])

    return SubspaceResult(
        q_new=q_new,
        harmonized_modes=modes_harm,
        singular_values=sigma_clip,
        principal_angles_deg=angles,
        omega=omega,
        rotation=r,
        matching_coefficients=r_m.copy(),
    )
