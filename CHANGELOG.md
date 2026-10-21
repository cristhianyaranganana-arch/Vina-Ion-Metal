# Changelog

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
