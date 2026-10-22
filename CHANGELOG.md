# Changelog

## Vina-Ion-Metal 1.0.1 (2026-10-23)

Maintenance patch: corrects five defects found by the independent
post-release audit of the metal coordination scoring path. The
byte-identical guarantee versus AutoDock Vina 1.2.7 (receptors without
parametrized metals) is preserved and remains CI-enforced.

### Correctness fixes (metal scoring path)

- **Coordination shell continuity (F1).** The coordination shell now
  extends to `cn_cutoff + 6*sigma` so the complete sigmoid tail of the CN
  weight stays inside the layer. Cutting exactly at `cn_cutoff` produced a
  discontinuity (~0.016 kcal/mol) when a donor crossed the shell border,
  contradicting the documented C-infinity smoothness of the CN penalty.
- **Metal pair force capping (F2).** The pair derivative is now smoothly
  saturated at 1e4 kcal/mol/A (`tanh`, the same scale as the standard
  `curl(v=1000)` cap). Near r -> 0 the derivative previously reached
  1e8-1e17 and could destabilize the BFGS local optimizer. Energies were
  already capped (AD4 convention); scores are unaffected.
- **`index_for_ad_type_with_state` contract (F3).** Now returns the exact
  entry for the requested oxidation state, or `ions.size()` when that state
  has no parametrized entry, as documented in the header. The previous
  silent first-match fallback made the unparametrized-state diagnostic of
  `Vina::check_receptor_metals()` unreachable (dead code). The scoring
  path is unchanged: it uses `index_for_ad_type()`, which keeps its
  fallback.
- **Dead cache field removed (F4).** `model_cache::key` was stored but
  never compared; cache validity is decided by the atom counts
  (`n_grid_atoms`/`n_movable`).
- **`non_cache` force parity (F5).** G0-G3 pseudoatom forces are now zeroed
  in `non_cache::eval_deriv`, matching `cache::eval_deriv`.

### Validation

- Metal gradient checks 4/4 PASS (analytic vs finite differences, both CN
  regimes and the no-metal control).
- CN detection, unparametrized-metal regression and the byte-identical
  regression against the official AutoDock Vina 1.2.7 binary all PASS
  (1iep: -13.170 both).
- Fixed-pose scores at the bundled examples are unchanged (glft2: 2.844).

### Compatibility

- Identical results to 1.0.0 for the validated poses. Docking trajectories
  can differ from 1.0.0 only in the rare configurations where a donor
  crosses the old shell border (bounded energy shift ~0.016 kcal/mol).

---

## Vina-Ion-Metal 1.0.0 (2026-10-21)

Corrected and completed successor of VinaMetal 1.0.0. Incorporates the fixes
for the 16 problems of the October 2026 technical audit plus additional
defects found during this revision. Based on AutoDock Vina 1.2.7.

### Critical fixes (P1)

- **Grid over-exclusion bug (audit #1).** `cache.cpp` excluded from the
  donor-type affinity maps every atom whose element mapped to the generic
  metal element (`EL_TYPE_Met`) — not only the parametrized ions. Receptors
  containing unparametrized metals (Na⁺, K⁺, Co, Cd, Hg, …) silently lost the
  standard Vina interaction between that metal and every ligand donor,
  breaking the byte-identical promise. Now only parametrized ions are excluded
  (`metal_scoring::is_metal_ion`); unparametrized metals keep the exact stock
  Vina 1.2.7 map treatment. New regression test
  `tests/test_unparametrized_metal.py` guards this behavior in CTest.
- **Scoring weights (audit #9).** The placeholder defaults (1.0, 1.0, 1.0)
  are replaced by provisional literature-informed heuristics: vdw 1.2,
  cn 0.5, geom 1.0 (rationale in `docs/USAGE.md`). A runtime notice is printed
  whenever coordination scoring is active. The full calibration benchmark
  remains pending and is documented in `LIMITATIONS.md` (unchanged
  recommendation: calibrate before virtual screening).

### Usability fixes (P2)

- **Dead option `--metal_cn_strict` implemented (audit #2).** The CN penalty
  now has two regimes: *lenient* (default) penalizes only over-coordination
  through a C-infinity softplus hinge — an incomplete ligand shell is not
  punished — and *strict* (`--metal_cn_strict`) penalizes any deviation from
  the target with the pure quadratic. Gradients implemented and covered by a
  new strict-mode gradient check. NOTE: VinaMetal 1.0.0 always behaved as
  strict (the flag was dead), so re-docking results may shift.
- **Cu²⁺/Ni²⁺ preparation path (audit #5).** `scripts/prepare_receptor.py` is
  now a self-contained pipeline that types parametrized ions directly into the
  PDBQT — including Cu/Ni, which Meeko cannot emit. Meeko (when installed)
  handles the protein part; the typed metal records are appended.
- **Receptor preparation robustness (audit #11).** The new pipeline resolves
  alternate locations, checks chain integrity (numbering gaps, missing
  backbone atoms), handles multi-model files, removes waters/HETATM on
  request, disambiguates metal records (CA ion vs alpha carbon) by residue
  identity, and prints a coordination pre-flight report. A dependency-free
  `--backend simple` converter (united-atom amino-acid templates) is included
  for the vina/vinardo scoring functions.
- **Occupancy detection hardened (audit #10).** Besides the typed donor count
  used for the CN target, a typing-independent geometric count of N/O/S atoms
  inside the coordination cutoff is computed; divergences are reported as
  preparation-quality diagnostics (both by the runtime
  `report_cn_targets()` and by the preparation script).
- **`--metal_oxidation` strict validation (audit #6).** Non-numeric values
  (`Fe=abc`) are hard errors instead of silently selecting the +2 entry;
  unknown elements and out-of-range states are rejected.
- **Oxidation-state fallback warning (audit #7).** When a declared state has
  no parametrized entry (e.g. `Fe=4`), a warning listing the valid states is
  printed at receptor load time.

### Project polish (P3)

- **Branding (audit #4).** Version string, issue-tracker URL and citation
  banner are Vina-Ion-Metal (`Vina-Ion-Metal v1.0.0 (based on AutoDock Vina
  1.2.7)`); log prefixes updated accordingly.
- **Dead variable `dE_dcos` removed (audit #3).**
- **Orphan `data/AD4Zn.dat` removed (audit #8).** It was a mislabeled copy of
  `AD4.1_bound.dat`, not the real AD4Zn parameter set; TODOs now point to the
  upstream AutoDock repository.
- **Portable binaries (audit #12).** `Boost_USE_STATIC_LIBS` is ON by default
  for source builds, with a graceful fallback (warning + dynamic linking)
  when static libraries are unavailable.
- **Regression test in CTest (audit #13).** `test_regression_no_metal.sh` is
  registered in CTest, locates the reference Vina 1.2.7 binary through
  `VINA_REFERENCE_BIN` (or skips with an explanatory message when absent),
  and the CI downloads the official reference binary so the byte-identical
  check runs on every push.
- **CI workflow fixed and hardened (audit #14 + extra finding).** The
  VinaMetal workflow contained a corrupted trigger (`branches: ain]`) that
  silently disabled CI on push/PR — fixed to `branches: [main]`. The Linux
  x86-64 scope is kept (intentional for v1.0); artifacts renamed to
  Vina-Ion-Metal.
- **AutoDock4Zn citation corrected (audit #17, extra finding).** The README
  cited *J. Chem. Inf. Model.* **2014**, 54, 1931–1941; the correct article is
  54(8), 2371–2379 (DOI 10.1021/ci500209e), with authors
  Santos-Martins, Forli, Ramos & Olson.
- **Recommended protocol documented (audits #15/#16).** `docs/USAGE.md` now
  spells out the multi-seed (≥3 seeds, exhaustiveness ≥16) and
  symmetry-aware-RMSD (RDKit `GetBestRMS`) protocol, plus PAINS counter-screen
  guidance for metal-chelating chemotypes.

### Documentation

- Bilingual README (English + `README.es.md`), updated `docs/USAGE.md`,
  `docs/LIMITATIONS.md`, `NOTICE` (attribution chain), `CITATION.cff`,
  `.zenodo.json` and packaging metadata for the Vina-Ion-Metal identity.

### Compatibility

- Receptors and ligands prepared for VinaMetal 1.0.0 work unchanged.
- Scoring without parametrized metals remains byte-identical to AutoDock
  Vina 1.2.7 (`--local_only` regression, CI-enforced).
- Scoring **with** parametrized metals changes with respect to VinaMetal
  1.0.0 (lenient CN default + heuristic weights + slightly stronger vdW
  term); re-run your experiments when comparing versions.

---

## VinaMetal 1.0.0 (2026-10-19)

Initial release. Fork of AutoDock Vina 1.2.7 with native metal coordination
scoring (AD4 12-6 metal pair potential, CN penalty, bounded geometry bonus,
analytic gradients) and automatic receptor occupancy detection. See
https://github.com/cristhianyaranganana-arch/VinaMetal and the October 2026
technical audit for the 16 issues corrected in Vina-Ion-Metal 1.0.0.
