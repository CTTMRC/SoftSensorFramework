# Release scope — 2026-10-09

This release extracts independently runnable functionality from the wider
research workspace. It provides PLS1, Groupwise PLS, elastic net, bagged trees,
numeric data fitting, compact inference, artificial examples, and opt-in saving.
The MATLAB namespace avoids collisions with historical helpers. MIT was selected
for redistribution; copyright is attributed to project contributors because the
available source did not establish a single legal holder.

The release excludes research data, private imports, original fitted models,
workbooks, study figures, result directories, historical archives, and local
review/handoff files. No source is imported from the confidential directory.
Historical factorial orchestration, its preprocessing/time-alignment combinations,
uncertainty recovery, and unreadable/incomplete legacy model entry points are
excluded. This package is not a drop-in replacement for historical scripts or a
reproduction of scientific results. Original research files remain outside it.

Changes in the extracted interfaces include fitting preprocessing within folds,
holding out whole batches, validating group and row alignment, using compact
saved inference models, reusing PLS coefficient paths where possible, and staging
unique output directories. No data or fitted model binaries ship in the archive.
Synthetic exports default to MAT; CSV/both require an explicit choice.

During this release's verification, an invalid seed-bound expression was repaired
in three argument declarations. Split masks now reject values other than zero/one;
unknown saved model families are rejected even for empty inference. Nonfinite
component/block settings are rejected. Stale generator/export help was corrected.

Tests and static checks establish the checked behavior on the recorded MATLAB
version. They do not establish absolute absence of bugs, scientific validity,
privacy of future caller outputs, or compatibility with every MATLAB version.
