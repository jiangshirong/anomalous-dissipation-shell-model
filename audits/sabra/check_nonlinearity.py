"""Audit the quadratic-energy balance of the repository's Sabra term.

This is intentionally a standalone audit.  It does not modify the cloned
repository and it does not silently replace the repository formula.

The audited expression is the one in
``functional-stability-theory/scripts/turbulence/compute_goy_shell_dfc.py``:

    N_n(u) = i * [ k_n     conj(u[n+1]) u[n+2]
                 - 1/2 k[n-1] conj(u[n+1]) u[n-1]
                 - 1/4 k[n-2] u[n-1] u[n-2] ].

For each state, the program computes

    D_NL = 2 Re sum_n conj(u[n]) N_n(u),

plus the contribution of every shell.  The isolated-triplet test activates
only three adjacent shells, so it is not explained away by a distant
boundary or by a long random state.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


REPO_FORMULA = (
    "i*(k_n*conj(u[n+1])*u[n+2] "
    "- 0.5*k[n-1]*conj(u[n+1])*u[n-1] "
    "- 0.25*k[n-2]*u[n-1]*u[n-2])"
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "sabra_audit_results"
DEFAULT_N_VALUES = (8, 12, 16)
N_RANDOM_STATES = 32
N_TRIPLET_STATES = 16
LAMBDA = 2.0
K0 = 1.0
TOLERANCE = 1.0e-12


def repo_nonlinear(u: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Exact vectorized transcription of the repository implementation."""
    n_shells = len(u)
    nl = np.zeros(n_shells, dtype=np.complex128)
    nl[: n_shells - 2] += k[: n_shells - 2] * np.conj(u[1 : n_shells - 1]) * u[2:n_shells]
    nl[1 : n_shells - 1] += (
        -0.5 * k[: n_shells - 2] * np.conj(u[2:n_shells]) * u[: n_shells - 2]
    )
    nl[2:n_shells] += -0.25 * k[: n_shells - 2] * u[1 : n_shells - 1] * u[: n_shells - 2]
    return 1j * nl


def shell_energy_rate(u: np.ndarray, nl: np.ndarray) -> np.ndarray:
    """Per-shell contribution to d(sum |u|^2)/dt from the nonlinear term."""
    return 2.0 * np.real(np.conj(u) * nl)


def relative_defect(total_defect: float, u: np.ndarray, nl: np.ndarray) -> float:
    scale = 2.0 * float(np.sum(np.abs(u) * np.abs(nl)))
    return abs(total_defect) / max(scale, np.finfo(float).tiny)


def random_state(rng: np.random.Generator, n_shells: int) -> np.ndarray:
    values = rng.normal(size=n_shells) + 1j * rng.normal(size=n_shells)
    values /= np.linalg.norm(values)
    return values


@dataclass
class RandomResult:
    n_shells: int
    seed: int
    nonlinear_energy_rate: float
    relative_defect: float
    interior_energy_rate: float
    boundary_energy_rate: float
    max_shell_abs_rate: float
    conserving_within_tolerance: bool


@dataclass
class TripletResult:
    n_shells: int
    triad_start: int
    seed: int
    nonlinear_energy_rate: float
    relative_defect: float
    shell_rates: str
    conserving_within_tolerance: bool


def audit_random_states() -> tuple[list[RandomResult], list[dict[str, object]]]:
    results: list[RandomResult] = []
    shell_rows: list[dict[str, object]] = []
    for n_shells in DEFAULT_N_VALUES:
        k = K0 * LAMBDA ** np.arange(n_shells, dtype=float)
        rng = np.random.default_rng(1000 + n_shells)
        shell_rates_all = []
        for seed in range(N_RANDOM_STATES):
            u = random_state(rng, n_shells)
            nl = repo_nonlinear(u, k)
            rates = shell_energy_rate(u, nl)
            total = float(np.sum(rates))
            boundary = float(rates[:2].sum() + rates[-2:].sum())
            interior = float(total - boundary)
            result = RandomResult(
                n_shells=n_shells,
                seed=seed,
                nonlinear_energy_rate=total,
                relative_defect=relative_defect(total, u, nl),
                interior_energy_rate=interior,
                boundary_energy_rate=boundary,
                max_shell_abs_rate=float(np.max(np.abs(rates))),
                conserving_within_tolerance=abs(total) <= TOLERANCE,
            )
            results.append(result)
            shell_rates_all.append(rates)

        shell_rates = np.asarray(shell_rates_all)
        for shell, values in enumerate(shell_rates.T):
            shell_rows.append(
                {
                    "n_shells": n_shells,
                    "shell": shell,
                    "mean_signed_rate": float(np.mean(values)),
                    "mean_abs_rate": float(np.mean(np.abs(values))),
                    "max_abs_rate": float(np.max(np.abs(values))),
                }
            )
    return results, shell_rows


def audit_isolated_triplets() -> list[TripletResult]:
    results: list[TripletResult] = []
    for n_shells in DEFAULT_N_VALUES:
        k = K0 * LAMBDA ** np.arange(n_shells, dtype=float)
        rng = np.random.default_rng(2000 + n_shells)
        for triad_start in range(1, n_shells - 3):
            for seed in range(N_TRIPLET_STATES):
                u = np.zeros(n_shells, dtype=np.complex128)
                local = random_state(rng, 3)
                u[triad_start : triad_start + 3] = local
                nl = repo_nonlinear(u, k)
                rates = shell_energy_rate(u, nl)
                total = float(np.sum(rates))
                results.append(
                    TripletResult(
                        n_shells=n_shells,
                        triad_start=triad_start,
                        seed=seed,
                        nonlinear_energy_rate=total,
                        relative_defect=relative_defect(total, u, nl),
                        shell_rates=json.dumps([float(x) for x in rates]),
                        conserving_within_tolerance=abs(total) <= TOLERANCE,
                    )
                )
    return results


def deterministic_counterexample() -> dict[str, object]:
    """A small exact-looking counterexample useful for manual inspection."""
    n_shells = 8
    k = K0 * LAMBDA ** np.arange(n_shells, dtype=float)
    u = np.zeros(n_shells, dtype=np.complex128)
    u[1:4] = np.asarray([1.0 + 0.0j, 1.0 + 0.0j, 0.0 + 1.0j])
    nl = repo_nonlinear(u, k)
    rates = shell_energy_rate(u, nl)
    total = float(np.sum(rates))
    expected = -3.5 * k[1]
    return {
        "n_shells": n_shells,
        "active_shells": [1, 2, 3],
        "u_active": ["1", "1", "i"],
        "k_active_start": float(k[1]),
        "per_shell_energy_rates": [float(x) for x in rates],
        "nonlinear_energy_rate": total,
        "-7/2*k_start": float(expected),
        "matches_-7/2*k_start": bool(np.isclose(total, expected, rtol=0.0, atol=1.0e-14)),
    }


def write_csv(path: Path, rows: list[object]) -> None:
    dictionaries = [asdict(row) if not isinstance(row, dict) else row for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def summarize(
    random_results: list[RandomResult],
    triplet_results: list[TripletResult],
    deterministic: dict[str, object],
) -> dict[str, object]:
    random_by_n: dict[str, object] = {}
    for n_shells in DEFAULT_N_VALUES:
        values = [r for r in random_results if r.n_shells == n_shells]
        defects = np.asarray([r.nonlinear_energy_rate for r in values])
        relative = np.asarray([r.relative_defect for r in values])
        random_by_n[str(n_shells)] = {
            "samples": len(values),
            "max_abs_nonlinear_energy_rate": float(np.max(np.abs(defects))),
            "mean_abs_nonlinear_energy_rate": float(np.mean(np.abs(defects))),
            "max_relative_defect": float(np.max(relative)),
            "all_within_tolerance": bool(np.all(np.abs(defects) <= TOLERANCE)),
        }

    triad_defects = np.asarray([r.nonlinear_energy_rate for r in triplet_results])
    summary = {
        "formula": REPO_FORMULA,
        "lambda": LAMBDA,
        "k0": K0,
        "tolerance": TOLERANCE,
        "random_state_summary": random_by_n,
        "isolated_triplet": {
            "samples": len(triplet_results),
            "max_abs_nonlinear_energy_rate": float(np.max(np.abs(triad_defects))),
            "mean_abs_nonlinear_energy_rate": float(np.mean(np.abs(triad_defects))),
            "all_within_tolerance": bool(np.all(np.abs(triad_defects) <= TOLERANCE)),
        },
        "deterministic_counterexample": deterministic,
        "interpretation_guardrail": (
            "A nonzero defect identifies failure of the audited expression for the "
            "stated finite-shell truncation; it does not by itself identify the "
            "intended corrected Sabra convention."
        ),
    }
    return summary


def write_markdown(summary: dict[str, object]) -> None:
    lines = [
        "# Sabra nonlinear energy audit",
        "",
        "This standalone audit reproduces the nonlinear expression in the downloaded repository.",
        "It does not modify repository code and does not assume a corrected convention.",
        "",
        "## Quantity tested",
        "",
        "`D_NL = 2 Re sum_n conj(u_n) N_n(u)`",
        "",
        f"Audited expression: `{REPO_FORMULA}`",
        "",
        "## Random-state results",
        "",
        "| shells | samples | max abs(D_NL) | mean abs(D_NL) | max relative defect | all within tolerance |",
        "|---:|---:|---:|---:|---:|:---:|",
    ]
    random_summary = summary["random_state_summary"]
    assert isinstance(random_summary, dict)
    for n_shells, item in random_summary.items():
        assert isinstance(item, dict)
        lines.append(
            f"| {n_shells} | {item['samples']} | {item['max_abs_nonlinear_energy_rate']:.6e} | "
            f"{item['mean_abs_nonlinear_energy_rate']:.6e} | {item['max_relative_defect']:.6e} | "
            f"{item['all_within_tolerance']} |"
        )

    triplet = summary["isolated_triplet"]
    assert isinstance(triplet, dict)
    counterexample = summary["deterministic_counterexample"]
    assert isinstance(counterexample, dict)
    lines.extend(
        [
            "",
            "## Isolated three-shell test",
            "",
            f"Samples: `{triplet['samples']}`",
            f"Maximum absolute defect: `{triplet['max_abs_nonlinear_energy_rate']:.6e}`",
            f"Mean absolute defect: `{triplet['mean_abs_nonlinear_energy_rate']:.6e}`",
            f"All within tolerance: `{triplet['all_within_tolerance']}`",
            "",
            "## Deterministic counterexample",
            "",
            "For `N=8`, `k_0=1`, with only `u_1=1`, `u_2=1`, `u_3=i` nonzero:",
            "",
            f"Per-shell rates: `{counterexample['per_shell_energy_rates']}`",
            f"Total `D_NL`: `{counterexample['nonlinear_energy_rate']:.15g}`",
            f"`-7/2*k_1`: `{counterexample['-7/2*k_start']:.15g}`",
            f"Exact numerical match: `{counterexample['matches_-7/2*k_start']}`",
            "",
            "## Interpretation boundary",
            "",
            summary["interpretation_guardrail"],
            "",
            "The per-shell CSV should be inspected together with the boundary columns.",
            "A later comparison against a separately specified standard Sabra convention is required before proposing a code correction.",
            "",
        ]
    )
    (OUT / "SABRA_NONLINEAR_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    random_results, shell_rows = audit_random_states()
    triplet_results = audit_isolated_triplets()
    deterministic = deterministic_counterexample()
    summary = summarize(random_results, triplet_results, deterministic)

    write_csv(OUT / "random_state_defects.csv", random_results)
    write_csv(OUT / "per_shell_defects.csv", shell_rows)
    write_csv(OUT / "isolated_triplet_defects.csv", triplet_results)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_markdown(summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"\nWrote audit outputs to: {OUT}")


if __name__ == "__main__":
    main()
