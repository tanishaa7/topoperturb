# TopoPerturb research status

**1 November 2026. The repository now contains the recovered earlier analyses and
preprocessed expression inputs. This closes most code-availability gaps, but does
not validate biological topology changes or establish exact historical replay.**

## What exists

- [src/](../src/) contains the canonical benchmark: input preparation, methods,
  synthetic/real comparisons and plotting. [src/earlier/](../src/earlier/) contains
  ten recovered modules for pilot PCA/Mapper, initial PH, Phase 1 robustness and
  interpretation, Phase 2 baselines/spectral PH, RPE1 and surrogates. The
  [README](../README.md) maps scripts to outputs and commands.
- [data/](../data/) has 13,534 K562 and 13,899 RPE1 cells in 30 PCs; metadata,
  counts and pilot targets; and recovered expression arrays with gene order:
  K562 13,534×622 (pilot plus controls), RPE1 91,891×383 (all non-excluded cells).
  [Provenance](../data/provenance.json) records hashes. These are preprocessed
  expression caches, not raw sequencing/count data or authenticated source revisions.
- [results/](../results/) holds 480 synthetic and 225 real comparisons, their two
  summary tables and benchmark figure. [archive.zip](archive.zip) retains earlier
  outputs, notebooks, reports and exact original recovered source versions.
  [requirements.txt](../requirements.txt) pins the verification environment.

## Pipeline and verified state

`prepare_data.py` filters `excluded`, selects targets and fits separate randomized
30-PC spaces per line, seed 42. It adds no normalization, log transform, gene
scaling or QC. K562 uses all retained pilot/control cells; RPE1 adds the strongest
eligible full-expression centroid shift as a selected pipeline positive control.
The pilot instead uses balanced groups. Cached loading now validates shapes,
genes, finite values and requested labels, and preserves array layout.

`synthetic.py` evaluates 12 geometry/ring/mixture conditions × 40 replicates;
`real.py` evaluates eight K562 targets, four shared RPE1 targets, one selected RPE1
control and two control-null groups × 15 draws. Each comparison uses 150 versus
150 cells. Draws sample without replacement internally but reuse cells across draws.

`methods.py` computes energy/covariance baselines, raw H1, normalized H0/H1,
corrected-resistance H1, diffusion H1 and three Gaussian-surrogate comparisons.
Rips PH uses F₂ and finite H0 bars. Separate symmetric kNN graphs start at 15 and
increase until connected; diffusion time is 8. Seven statistics use 199 pooled
label permutations with plus-one p-values; three use 99 Gaussian surrogate pairs.
Seed tuples are `[42, condition_index, replicate]` and `[42,77,group_index,draw]`.
BH applies to primary real draw 0 across genes per method and cell line.

**Executed evidence** ([verification.json](verification.json)):

- 14 tests passed; two expected failures document small-cloud kNN and zero-scale
  handling. All current table aggregates/BH values match saved draws. Previous
  five-comparison benchmark replays and figure reconstruction matched the archives.
- Recovered expression arrays were verified value-for-value after compression.
  Rebuilding PCA selects identical cells but differs numerically by about 0.002
  at most in coordinates; original PCA clouds were preserved.
- The recovered pilot reproduces saved descriptive statistics within about 2×10⁻⁶.
  Selected original 2,000-permutation Phase 1, normalized PH, conventional and RPE1
  checks reproduce their p-values. Some statistics exceed the strict replay
  tolerance (maximum checked discrepancy about 5×10⁻⁴); the checker reports failure
  rather than hiding this. Original environment/platform provenance is incomplete.
- Full earlier grids and the initial hash-seeded analysis have not been rerun.
  Generators for the old scale-test summary and combined detection-matrix
  assembly remain missing; the source analyses feeding them are now present.

## What the results support

[Synthetic results](../results/synthetic.csv): raw H1 rejects scaling on 40/40 and
anisotropy on 27/40 replicates; normalized H1 rejects on 10/40 and 30/40. For the
K562-spectrum ellipse, raw/normalized/diffusion H1 reject on 1/40, 4/40 and 27/40.
These demonstrate metric-dependent diagram sensitivity, not direct loop recovery.

[Real results](../results/real.csv): primary diffusion q=0.04 for K562 RPL21 and
q=0.02 for RPE1 ZMAT2/CPSF3; rejection fractions are 13/15, 13/15 and 15/15.
Control-null diffusion rejects on 2/15 K562 and 0/15 RPE1 draws; surrogate H0
rejects on 5/15 K562 controls. All eight K562 targets have energy rejection on
every draw. Additional biological information from PH has not been established.

## Remaining issues

1. **Inference:** pooled permutations test distribution equality, not topology
   modulo geometry. Normalized statistics do not guarantee a normalized null.
   Low rejection under one anisotropy model does not prove effects beyond covariance.
   Gaussian plug-in comparisons and non-rejection do not establish equivalence.
2. **Synthetic truth:** noisy rings and Gaussian mixtures generally have full
   support; they model latent/density structure. Covariance matching is in
   expectation. The figure's “known topology / should stay ~5” wording needs revision.
3. **Calibration/replication:** minimum benchmark p is 0.005; RPL21's q uses that
   floor. BH does not cover all methods/exploratory choices. Shared cells and a
   same-data selected positive control do not provide independent replication.
   No guide/batch metadata or biological loop validation is available.
4. **Recovered implementation:** initial PH seeds depend on randomized Python
   `hash()`; its original hash seed is unknown. Pilot nulls are 75-versus-75 even
   though the original output labelled pilot size 150; the adapted output makes
   this explicit. Centered/scaled MMD fixes a bandwidth chosen after label-dependent
   transformations, weakening exact permutation calibration. Representative cycles
   are exploratory reconstructions, not certified biological loops.
5. **Reproducibility/reporting:** input source revisions, normalization and symbol
   annotation provenance remain unknown. Some PCA/statistic replays differ slightly;
   matrix-arithmetic warnings are recorded. Small/degenerate-cloud edge cases remain.
   Plot titles still hard-code 40/15. Archived claims that dispersion detects all
   genes or that every robustness setting is BH-significant overstate their tables.

## Contribution and next steps

The contribution is an empirical perturbation-comparison benchmark of geometric
confounding and metric-dependent sensitivity. Geometry sensitivity and spectral
PH are prior work; no new PH algorithm, validated biological topology detector or
“first” claim is established.

**Essential for the poster, in order**

1. [ ] Use that narrower contribution; correct figure truth labels and covariance
   claims; report denominators, uncertainty, correction scope and both control rates.
2. [ ] Freeze a repository revision and inputs; document recovered-data provenance
   and numerical replay differences. Use the saved PCA clouds for the reported benchmark.
3. [ ] Before strengthening real-data claims, run a bounded, predeclared check of
   candidate/control permutation resolution and PCA/k/time sensitivity.
4. [ ] If claiming loop detection, validate persistent features against matched
   ring/mixture controls and biological metadata; otherwise report diagram differences.

**For a benchmark paper**

1. [ ] Define the target—distribution change, latent-loop recovery or topology up
   to geometry—and establish appropriate nuisance-aware calibration.
2. [ ] Separate selection from evaluation; broaden effect/noise settings, seeds,
   uncertainty and explicit feature-recovery benchmarks.
3. [ ] Obtain guide/batch information and independent replication; demonstrate
   interpretable value beyond conventional distributional baselines.
4. [ ] Resolve historical replay/provenance, complete summary assembly and full
   reruns, repair calibration weaknesses and review related work before generality claims.
