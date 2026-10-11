# Verification record

Date: 2026-10-09 (Europe/Rome). Runtime: MATLAB R2025a Update 1 on Windows,
with Statistics and Machine Learning Toolbox.

The current release passed all 67 class-based tests, with no failed or incomplete
tests. Code Analyzer reported no messages across 13 public functions and three
test classes. The suite checks:

- Exact/OLS/independent `plsregress` agreement, scores, rank and constant cases.
- CV coefficient-prefix agreement with separate fits, predictor-index expansion,
  and validation-label mutation with scaling, block balancing, and reduction.
- Both artificial targets across all four families and repeatable caller RNG.
- Batch disjointness, arbitrary/unseen groups, paired missing-row removal,
  invalid masks/shapes, empty inference, and metric edge cases.
- MAT save/reload predictions for all four families, including compact trees.
- No default persistent writes, unique saves preserving prior files, MAT/CSV/both
  choices, and write failures leaving only uncommitted `.partial` directories.

The freshly extracted archive passed all 67 tests in a separate MATLAB process
with only MATLAB and the release root on the path. It had zero dependencies on
the surrounding research workspace. An additional restricted run passed all 52
base tests using only MATLAB's base/shared paths, with `lasso` and `plsregress`
unavailable. All 67 tests also passed with only base/shared/Statistics paths.
Thus System Identification Toolbox is not a runtime requirement; the static
scanner's extra product match came from the globally visible `predict` name.

The package uses an explicit text-file allowlist and SHA-256 manifest. The local
audit checked 645 excluded assets by size and, where needed, exact hash; it found
no exact asset matches, original asset filenames, forbidden paths/imports, or
selected confidential workbook identifier matches. All distributed files are
UTF-8 text. There are no data/model binaries, workbooks, figures, result exports,
or private review artifacts in the archive. These bounded checks do not prove
the absence of every conceivable private value. Research arrays were not loaded.

[VERIFICATION.json](VERIFICATION.json) records test results, restricted-path
checks, audit counts, and hashes of the verified MATLAB source/test files. Final
package rebuilding changes only documentation/evidence; the MATLAB files remain
byte-identical to those tested in isolation. `MANIFEST.sha256` covers all other
distributed files; the archive has a separate SHA-256 sidecar.

No original research arrays or study-result claims are used by these tests.
Other MATLAB versions and operating systems were not tested. This record does
not establish absolute bug-freedom or reproduction of the scientific study.
