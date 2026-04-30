# TikZ Train-Track Plotter Plan

## Goal

Create a first reliable plotting pipeline from train-track coding to vector output, with prong-level tangency represented explicitly, and with direct PDF generation through LaTeX.

## Phase 0 - Scope and Constraints

- Build a plotting tool inside the existing C++/CMake workflow.
- Prioritize combinatorial correctness and tangency representation over global layout optimization.
- Keep output deterministic and reproducible for paper figures.
- Avoid introducing new mandatory external runtime dependencies.

## Phase 1 - First Deliverable (TikZ-first)

- Add `examples/ttplot.cpp` as a command-line program.
- Accept coding input from:
  - `--coding "..."` for inline coding,
  - `--coding-file <path>` for file-based input,
  - stdin fallback when no coding flag is provided.
- Parse both coding forms used in this project:
  - modern labeled block form (5 digits per block, 1-based print convention),
  - legacy unlabeled block form (4 digits per block), mapped to label `0` internally.
- Reconstruct the track using `traintrack(const intVec&)`.
- Validate topology with `check()`.

## Phase 2 - Geometry and Tangency Model

- Place multigon centers on a deterministic circle layout.
- For each multigon:
  - assign prong directions,
  - assign one endpoint slot per edge attached to each prong,
  - offset endpoint slots tangentially to represent cusp/tangency structure.
- Draw each edge as a cubic Bezier between endpoint slots.
- Use endpoint control points aligned with prong direction to enforce local tangency.

## Phase 3 - Output and Labeling

- Emit TikZ in one of two modes:
  - standalone LaTeX document (`.tex`) for direct `pdflatex` to PDF,
  - snippet mode (`--snippet`) for direct inclusion in existing papers.
- Add label modes:
  - `none`, `multigons`, `prongs`, `edges`, `all`.
- Default style aimed at paper usage, with tunable scale/curvature parameters.

## Phase 4 - Validation

- Build target with CMake auto-discovery (`example_ttplot`).
- Run smoke tests on small codings.
- Verify deterministic output for repeated runs.
- If TeX is available, compile to PDF and visually inspect tangency at multi-edge prongs.

## Phase 5 - Documentation

- Add concise usage section in `README.md`.
- Document accepted coding formats and limits.
- Provide one or two example commands for quick figure generation.

## Near-Term Enhancements (v1.1-v1.3)

- Improve global layout:
  - crossing-reduction heuristics,
  - force-based refinement seeded from circular placement,
  - optional pinning of selected multigons.
- Improve labeling:
  - collision avoidance,
  - leader lines,
  - selective edge-label subsets.
- Add export helpers:
  - JSON geometry export,
  - SVG backend or converter path.

## Medium-Term Options

- Add Mathematica function in `mathematica/TrainTracks.m` that plots from coding directly.
- Add a Python plotting helper that uses the same geometry model and writes SVG/PDF.
- Expose a reusable geometry builder API in the C++ library for multi-backend plotting.

## Long-Term Option

- Create a backend-neutral plotting layer with interchangeable emitters (TikZ, Mathematica, Python/SVG), all driven by one canonical endpoint/tangency representation.

## Known Initial Limitations

- Compact block parsing is practical only when per-field values are single-digit.
- Initial layout is deterministic but not crossing-optimal for large tracks.
- Canonical edge ordering labels are approximated from coding traversal-derived weights in v1.
