# Changelog

All notable changes to `aoutools` are recorded here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the library is pre-1.0, breaking changes may land in a minor version.

## [Unreleased]

This release corrects how the effect allele is matched against the *All of Us*
VDS. **Scores computed with 0.1.2 or earlier can change** — for most files they
move because previously dropped variants now contribute; where the effect allele
was the reference base, earlier rankings could have been wrong, not merely
shifted. See **Fixed** below for who is affected.

### Added

- `init_hail`, `get_vds_path`, `get_workspace_bucket`, and `get_google_project` —
  Workbench bootstrap helpers that wire up requester-pays billing, the GRCh38
  reference, and the VDS/bucket/project paths in one place, so a notebook no
  longer pastes that setup by hand.

### Deprecated

- `read_prscs` — the name implied a tie to the PRS-CS tool, but it only applies
  one fixed header-less column layout. Call `read_prs_weights` directly with
  `header=False` and the PRS-CS `column_map`; `read_prscs` still works, warns, and
  will be removed in a future release.

### Removed

- **(Breaking)** The `PRSConfig` options `split_multi` and `ref_is_effect_allele`.
  Both selected a scoring path that silently lost data: `split_multi=False` zeroed
  every homozygous-reference sample at a reference-effect variant, and
  `ref_is_effect_allele` declared the effect allele's orientation for a whole file
  when it is really a per-variant property. Passing either now raises `TypeError`.

### Fixed

- Reference-effect weights are no longer silently dropped. Effect-allele
  orientation is now resolved per variant against the VDS reference, and variants
  are matched on the unordered allele set, so a SNP whose reference base sorts
  after its alternate, and an already-biallelic non-minimal variant, both match
  now instead of being skipped (922 of 1,940 variants were dropped for one real
  score before this fix). Scores that were computed from a biased subset of
  variants now include the full set.
- At a multi-allelic site, a reference-effect weight no longer over-credits a
  sample that carries a *different* alternate allele as though it were homozygous
  reference. This was per-sample and genotype-dependent, so it could reorder a
  cohort; anyone scoring reference-effect weights at multi-allelic sites is
  affected.

> **Note:** releases 0.1.0–0.1.2 targeted the *All of Us* Researcher Workbench
> 1.0, which was decommissioned on June 30, 2026. They are kept here for the
> record; current development targets Researcher Workbench 2.0 (new VDS path and
> Hail setup — see **Added** above).

## [0.1.2] - 2026-02-06

### Changed

- Raised the minimum supported Python version.

## [0.1.1] - 2025-10-11

- Maintenance release. (Predates this changelog; see the Git history for detail.)

## [0.1.0] - 2025-08-27

- Initial release: the `aoutools.prs` submodule — a flexible reader for PRS
  weight files and a cost-efficient strategy for calculating PRS directly on the
  *All of Us* VDS, including batch scoring.

[Unreleased]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1.2...dev
[0.1.2]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1...v0.1.1
[0.1.0]: https://github.com/dokyoonkimlab/aoutools/releases/tag/v0.1
