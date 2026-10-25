# Anomalous dissipation in shell models: two research results

[简体中文](README.zh-CN.md)

This repository presents two results. In a deterministic positive two-channel shell network, coefficient changes supported on a set of zero natural density can suppress anomalous dissipation. In a scalar shell chain with random coefficients, identical complete time laws at individual edges can coexist with opposite dissipation outcomes because dependence between edges differs. Both studies specify the equations, solution classes, initial and boundary conditions, and actual viscous loss, with complete analytical proofs.

## Result one: dissipation under zero-density modifications

The model is a forced positive two-channel network with zero initial data and viscous weights K_n², where K_n=2^n. Two configurations share the graph, forcing, initial data, viscosity and natural density of weak layers. The periodic exchange configuration has a viscosity-uniform positive loss bound on a sufficiently long fixed interval. Weakening exchange within consecutive blocks supported on a set of zero natural density makes loss vanish on every fixed interval. Relative coefficient changes at the modified layers can be large.

- [Model, main result and reading route](docs/sparse-modification/model-zh.md)
- [Complete proof of the main theorem](docs/sparse-modification/proofs-zh.md)
- [Further results: weak blocks, frozen arrangements, permutations, free decay and stability](docs/sparse-modification/extensions-zh.md)
- [Logarithmic exchange: proved results and remaining gap](docs/sparse-modification/logarithmic-zh.md)
- [Literature and applicability of existing theorems](docs/sparse-modification/literature-zh.md)

## Result two: sharp criteria in correlated random environments

The second model is an unforced positive scalar chain with a unit pulse at the root. Each edge switches between weak and strong coupling. Common resets and private flips inside independent finite blocks control dependence. Holding the complete stationary path law of each individual edge fixed still allows either zero or positive actual vanishing-viscosity loss.

For example, take blocks of j edges with single-edge correlation rate j². Independent edges give vanishing loss on every fixed interval; any fixed positive common-reset proportion gives anomalous dissipation. For arbitrary finite block lengths ℓ_j, positive finite correlation rates r_j and reset proportions θ_j∈[0,1], the classification uses the cost series

```math
\mathcal C_j=\begin{cases}
\ell_j/r_j,&\theta_j=0,\\
\min\{\ell_j,\theta_j^{-1}\}/r_j+(\mathcal N_j-1)/(\theta_jr_j),&0<\theta_j<1,\\
1/r_j,&\theta_j=1,
\end{cases}
\qquad \mathcal C=\sum_j\mathcal C_j.
```

Here 𝒩_j is the geometric iteration count defined by the model. Divergence gives almost-sure vanishing loss simultaneously on all finite intervals. Convergence gives positive expected loss on a deterministic interval and a common positive-probability pathwise lower bound. The almost-sure statement allows sample-dependent windows: taking the vanishing-viscosity lower or upper limit first, and then increasing time, recovers the full initial energy.

- [Model, environment and general classification theorem](docs/correlated-environments/model-zh.md)
- [Complete proof of the classification](docs/correlated-environments/proofs-zh.md)
- [Corollaries, geometric counterexamples and synchronization recovery](docs/correlated-environments/consequences-zh.md)
- [Supplementary methods and results beyond the general theorem](docs/correlated-environments/supplementary-methods-zh.md)
- [Literature comparison and the mathematical contribution](docs/correlated-environments/literature-zh.md)

Under this project's fixed criteria, D1 achieves the objective of a substantially original shell model with provable anomalous dissipation. The controlled comparisons, uniform estimates and comparisons with existing theorems are documented in the literature discussion. The scope and supporting evidence of that assessment remain explicit.

## Other proved results and open problems

- [Autonomous rotation chain: dissipation proof and rigorous three-shell backflow certificate](docs/other-results/rotation-chain-zh.md)
- [Autonomous phase gates: a fast/slow contrast at equal controller energy](docs/other-results/phase-gates-zh.md)
- [Earlier models and valid technical derivations](docs/other-results/earlier-models-zh.md)
- [Partial results, obstructions and stopping points](docs/other-results/partial-results-zh.md)
- [Current open problems](docs/open-problems-zh.md)

## Proofs, checks and sources

The mathematical conclusions rely on the analytical proofs in the research documents. Programs check finite identities, finite certificates and repository structure; their output states the scope actually checked. Literature applicability is tied to specific equations and theorems that were read. Technical documents are in Chinese; this README and its Chinese counterpart provide matching entry points.

```bash
python -m pip install -r requirements.txt
python tests/check_repository.py
python tests/check_identities.py
python tests/check_cut_obstruction.py
python tests/check_rotation_chain.py
python tests/check_repository_negative.py
```

The complete rotation certificate uses exact rational arithmetic and takes longer than the algebra checks; use --algebra-only for the latter. Default repository checks do not require local manuscripts or downloaded papers.

- [Reference index](references/README.md)
- [Key corrections and limits of the claims](docs/history/corrections-zh.md)
- [Research history](docs/history/research-history-zh.md)
- [Migration index for earlier sections](docs/history/migration-map-zh.md)
- [Organization and maintenance rules](docs/history/organization-zh.md)

## Local materials and licensing

Manuscripts, submission letters, editorial correspondence, third-party full texts and build tools live in the ignored local directory. The public repository includes the research results and proofs used in the papers, while manuscript files and submission materials stay local.

Original research documents are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt), and original programs under [MIT](LICENSES/MIT.txt); see the [license scope](LICENSE). Third-party materials retain their own licenses.
