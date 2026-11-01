# Methods: formulas, references, and how they are checked

Every method lives in a small module of the `topoperturb` package and is verified
against an independent reference by
[`experiments/verify_methods.py`](../experiments/verify_methods.py).

## Three questions, kept separate

A single "significant" result is easy to over-read. We separate three questions and
report them with different evidence:

| Question | Evidence | Where |
|---|---|---|
| Do the two clouds differ at all? | permutation test of any statistic | `twosample.permutation_test` |
| Does the difference exceed a Gaussian? | covariance-matched surrogate test | `twosample.surrogate_test` |
| *What* shape changed, and by how much? | leading H1 lifetime, signed change, bootstrap CI | `features.py` |

A significant diagram distance is a starting point, not proof of a biological loop.

## Spectral distances ([`distances.py`](../topoperturb/distances.py))

Persistent homology sees only the distance matrix it is given. On a kNN graph of
the cells, two intrinsic distances recover loops and clusters that Euclidean
distance blurs, following **Damrich, Berens & Kobak, NeurIPS 2024**:

- **Diffusion distance** at time `t`: distance between the landing distributions of
  `t`-step random walks, from the eigendecomposition of the normalised adjacency.
- **Effective-resistance distance**: resistance of the graph viewed as a network of
  unit resistors, with the von Luxburg correction that removes the degree term.

*Check:* the effective-resistance machinery (Laplacian pseudoinverse) matches
`networkx.resistance_distance` to ~1e-16; the correction is a documented additive
term on top of that.

## Persistent homology ([`homology.py`](../topoperturb/homology.py))

`ripser` builds the Vietoris–Rips persistence (H0 = components, H1 = loops) from a
distance matrix; each diagram point is a `(birth, death)` pair. Distances are scaled
by their mean off-diagonal entry so a pure size difference does not look like a
shape difference. Two diagrams are compared with the Wasserstein distance from
`persim`.

*Check:* a diagram's Wasserstein distance to itself is exactly zero.

## Two-sample tests ([`twosample.py`](../topoperturb/twosample.py))

- **Energy distance** (Székely & Rizzo): `2·E‖A−B‖ − E‖A−A'‖ − E‖B−B'‖`, a strong
  conventional comparator (also central to scPerturb). Zero only for identical
  distributions.
- **Permutation test**: shuffle the pooled cells `P_PERM = 199` times; the p-value is
  `(1 + #{shuffled ≥ observed}) / (1 + P_PERM)`. One shared set of shuffles tests
  every statistic. Raising the count improves *resolution*, not *validity*.
- **Surrogate test**: compare the observed statistic to the same statistic between
  Gaussian clouds matched to each group's mean and covariance (`M_SURR = 99`).
  Rejection means the difference is bigger than mean and covariance alone produce —
  it does **not** by itself prove a loop; any non-Gaussian structure can cause it.

*Check:* the energy distance matches its definition; the permutation test on
identical Gaussians rejects at about the nominal 5%.

## Interpretable features ([`features.py`](../topoperturb/features.py))

- **Leading H1 lifetime** — how persistent the main loop is ("is there a ring").
- **Signed change** — the leading lifetime in the perturbed cloud minus the control;
  its sign says whether a loop *appeared/strengthened* or *weakened/vanished*, which
  a distance alone cannot.
- **Bootstrap interval** — resample the cells to put a 95% band on the leading
  lifetime, so a reported loop always comes with its reproducibility.

## Calibration ([`experiments/04_calibration.py`](../experiments/04_calibration.py))

Because the representation and the normalised statistics are chosen before testing,
calibration is run as its own stage: comparing controls with controls across group
sizes, reporting each method's false-positive rate with a Wilson interval and a KS
test of p-value uniformity. Treat surrogate results as exploratory until a method
passes this check. Cells are drawn only from the testing controls, so the PCA axes
are independent of the null being measured (`docs/DATA.md`).

## Input space ([`experiments/05_compare_spaces.py`](../experiments/05_compare_spaces.py))

Nothing after step 1 needs PCA: `distances` and `homology` take any point cloud, so
the normalised gene matrix can be fed in directly. `05_compare_spaces` runs the same
battery on the same cells and permutations in 30-PC space and in gene space, and
reports graph connectivity, distance concentration, 15-NN agreement, null rate,
sensitivity and runtime. The 30-PC space is primary (better distance contrast, faster,
no null inflation); gene space is kept as a sensitivity baseline. Leading-H1 lifetimes
are representation-dependent, so a loop readout must always name its space. Results
and limits: `docs/UPDATE_REPORT.md` §4.
