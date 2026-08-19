# Changelog

All notable changes to `aoutools` are recorded here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the library is pre-1.0, breaking changes may land in a minor version.

## [Unreleased]

## [0.2.0] - 2026-08-19

This release resolves effect-allele orientation per variant, so a weights file no
longer has to be harmonized to a single orientation before scoring. **Scores
computed with 0.1.2 or earlier can change**: rows whose orientation did not match
the file-wide setting were skipped, so a score was computed from part of the file
rather than all of it, and sample rankings can shift accordingly. See **Removed**
and **Fixed** below for what this means for results you already have.

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

- **(Breaking)** The `PRSConfig` options `split_multi`, `ref_is_effect_allele`,
  and `strict_allele_match`. Passing any of them now raises `TypeError`.

  `ref_is_effect_allele` told the scorer that a whole file used the reference
  allele as its effect allele — it assumed you had harmonized the file to one
  orientation and were declaring which one. Nothing needs declaring now:
  orientation is resolved per variant against the VDS, so a file carrying both
  orientations, as PGS Catalog files commonly do, scores in full with no setting
  at all. If you used `ref_is_effect_allele=True` on a file that did meet the
  assumption, your **rankings were correct** — every score was short by the same
  constant, so percentile and z-score results still stand.

  `split_multi=False` selected an alternative scoring path that has been removed,
  and `strict_allele_match` tuned the allele check on that path only, so it went
  with it. Neither was the default.

### Fixed

- A weights file no longer has to be uniformly oriented. Effect-allele
  orientation is resolved per variant against the VDS reference, and variants are
  matched on the unordered allele set, so a SNP whose reference base sorts after
  its alternate, and an already-biallelic non-minimal variant, both match now
  instead of being skipped. Rows that disagreed with the file-wide setting used
  to be passed over quietly — not scored, and not counted in `n_matched` — which
  made a partial score hard to notice. Those rows now contribute.

- Alleles written in lowercase are now read correctly. Allele comparisons are
  case-sensitive, so a weights file using `a`/`g` rather than `A`/`G` did not
  match any variant. What you saw depended on a setting: with allele validation
  off — the default — the file loaded and every sample scored zero, with no
  error and no warning; with it on, as `calculate_pgs` always uses, the
  lowercase rows were dropped as invalid, and a wholly lowercase file was
  rejected outright. Alleles are now standardized to uppercase as a file is
  read, so all of these now score, and both the scores and `n_matched` change
  accordingly. PGS Catalog files are uppercase and are unaffected; this matters
  for weights files you assemble yourself.

  One consequence worth knowing: a file listing the same variant twice at one
  position, differing only in letter case (`a`/`g` and `A`/`G`), is now reported
  as a duplicate and rejected. Previously the lowercase row was quietly ignored
  and the file scored. The two rows genuinely conflict, so this is now an error
  you can see rather than a choice made for you.

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

[Unreleased]: https://github.com/dokyoonkimlab/aoutools/compare/v0.2.0...dev
[0.2.0]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/dokyoonkimlab/aoutools/compare/v0.1...v0.1.1
[0.1.0]: https://github.com/dokyoonkimlab/aoutools/releases/tag/v0.1
