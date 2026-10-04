# Anomalous dissipation in positive shell networks

**English** | [简体中文](README.zh-CN.md)

This project studies positive, conservative shell networks with scale-dependent exchange. Each scale has two nonnegative real amplitudes, w_n and z_n, with wave number K_n=2^n and energy (w_n²+z_n²)/2. The network is a mathematical cascade model. It is neither standard Sabra nor a rigorously derived reduction of three-dimensional Navier–Stokes.

The question is whether actual viscous energy loss remains positive as viscosity tends to zero on a fixed observation window. The infinite-shell solution is constructed at each positive viscosity before taking the vanishing-viscosity limit. Closed finite truncations conserve energy in the absence of viscosity and forcing; terminal leakage and numerical errors do not count as dissipation anomaly.

## Main result

With common zero initial data, unit forcing at both roots and viscosity weights K_n², the cumulative loss is

```math
Q_\nu(T)=\nu\int_0^T\sum_{n\ge0}K_n^2
\bigl((w_n^\nu)^2+(z_n^\nu)^2\bigr)\,dt,
\qquad E^\nu(T)+Q_\nu(T)=\int_0^T(w_0^\nu+z_0^\nu)\,dt.
```

A periodic exchange configuration has a positive, viscosity-independent lower bound for Q_ν(T)/T on sufficiently long fixed windows. Weakening coefficients only in sparse consecutive blocks suppresses the loss: Q_ν(T)→0 for every fixed T>0. The two configurations have the same graph, forcing, initial data and viscosity weights; both have weak-layer density 1/2, and their coefficients differ on a set of natural density zero. No edge is removed, and all exchange coefficients still grow with the shell index.

The modified coefficients undergo large relative changes. This result does not say that every zero-density modification suppresses dissipation. Statements and complete proofs are in [the research results, §7.6–7.8](docs/03-results-zh.md#76-密集弱交换逐层级数发散却仍有异常耗散); the actual internal review scope is recorded in [the evidence, §8.11](docs/04-evidence-zh.md#811-选定动态对照的证明审查与依赖闭合).

## Further results and open questions

| Topic | Current result |
|---|---|
| Short weak blocks | For the specified family N_j=2^(2^j), ℓ_j=⌈c log₂N_j⌉ with c>0, c<1/2 gives an anomaly and c≥1/2 gives vanishing fixed-window loss. |
| Frozen independent binary coefficients | With weak-layer probability p, p<1/4 gives an anomaly almost surely; p≥1/4 gives vanishing fixed-window loss almost surely. Coefficients are sampled once and then fixed; windows may depend on the realization. |
| Permuting actual exchange coefficients | Admissible permutations with displacement at most c log₂(n+2)+C preserve the anomaly for c<1/2. At the critical scale, a zero-density permutation with vanishing loss exists; not every critical permutation is classified. |
| Unforced decay | With common initial amplitudes 1 at both roots and zero elsewhere, the periodic configuration loses more than half the initial energy on one viscosity-independent window for all sufficiently small viscosities, while the selected critical permutation has vanishing loss on every fixed window. |
| Logarithmic exchange h_n=η2^(n/2)(n+1)^α | Fixed-window loss vanishes for α≤1/2; sufficiently long fixed windows have a positive lower bound for α>3/4. The interval 1/2<α≤3/4 remains open. |
| Loss by channel | In the logarithmic family with η=β=1 and common zero initial data, the z-channel viscous loss tends to zero for every fixed α>0. Any possible anomaly is carried by w. |

The general dynamic role of a divergent geometric-profile sum P remains open. In particular, α=2/3 is not classified: endpoint and receiver-occupancy conditions have not been established for the actual trajectory. The proposed fixed-viscosity continuity and zero-viscosity discontinuity task in [design §19](docs/05-design-zh.md#19-2026-10-04-选题建议固定黏性稳定与零黏耗散不连续) remains a proposal, not a proved theorem.

These are internal analytic results with explicitly recorded verification scopes. The conservative pump unit and several proof tools have known precedents. The concrete dynamic contrast is the contribution candidate; global priority, substantive originality and peer-review acceptance have not been certified. The core result has been submitted to JSP by the author, which does not establish acceptance.

## Reading guide

Start with [RESEARCH_HANDOFF.md](RESEARCH_HANDOFF.md) for the current state. The research documents are in Chinese.

| File | Purpose |
|---|---|
| [Overview](docs/01-overview-zh.md) | Variables, energy transfer and result boundaries |
| [Foundations](docs/02-foundations-zh.md) | Equations, solution class, construction and exact budgets |
| [Results](docs/03-results-zh.md) | Complete proofs and their assumptions |
| [Evidence](docs/04-evidence-zh.md) | Corrections, counterexamples, finite checks and literature comparisons |
| [Design](docs/05-design-zh.md) | Decisions, open tasks and stopping conditions |
| [Process](docs/06-process-zh.md) | Historical provenance and migration records; consult when needed |
| [Legacy appendix](docs/appendices/legacy-models-zh.md), [mechanisms appendix](docs/appendices/mechanisms-zh.md) | Earlier valid derivations and mechanism screening |
| [References](references/README.md) | Sources and the portions actually read |

## Repository and local files

GitHub contains the two READMEs, research handoff and conventions, research documents, reference index, finite-check programs and dependency lists.

The following stay local: `manuscript/`, `submission/`, `output/`, `tmp/`, `.local-backup/`, downloaded full papers, third-party checkouts and raw experimental data. This includes all manuscript sources, paper PDFs, cover letters and submission metadata. Research statements and proofs used in the manuscript may be synchronized in `docs/`.

This public repository starts with a reviewed research snapshot published on 2026-10-04. Earlier Git history, including the old W/S manuscript `manuscript/model.tex`, is preserved separately in a private history repository and a local backup; it is not part of this public repository. Earlier commit hashes in the process records refer to that preserved history. This publication does not imply that the work was never previously public. The synchronization allowlist and checker enforce the file boundary. A clone can read and check the research without access to the local manuscripts.

## Verification

Use Python 3.11 and the dependencies in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python tests/check_identities.py
python tests/check_cut_obstruction.py
python audits/logarithmic/check_review_claims.py
python tests/check_repository.py
python tests/check_repository_negative.py
```

Before an explicitly authorized upload, `python tests/check_repository.py --staged` checks the actual Git index. Authors who retain the excluded local material can additionally run `python tests/check_repository.py --local-manuscripts` for the existing bilingual manuscript checks; those files are not required by the default research check.

Finite algebra checks do not prove infinite-dimensional dynamics or originality. Repository checks validate file scope and navigation, not mathematical truth. There is no published complete time-integration solver interface. The separate [third-party Sabra audit](audits/sabra/README.md) is historical work and is not a contribution of this model.
