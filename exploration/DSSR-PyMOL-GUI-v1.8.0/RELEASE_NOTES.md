# 1.8.0 — sequence, letter colors and rectangle selection, 2026-09-06

One Python runtime file, one shared PyMOL Qt window. All scientific layout
implementations and command interfaces are retained.

## Requested interaction changes

- A/C/G/T/U have distinct letter colors, enabled by default. Unknown/other
  letters stay neutral. The existing neutral gel/flat node fills, white
  background and edge styles are unchanged. Options can turn letter colors off.
- A compact horizontal 1D sequence sits above the canvas, with original base
  symbols, chain separators and actual residue IDs. It shares colors and
  selection highlights with 2D/3D; clicking or dragging selects a range.
- Blank-canvas left-drag draws a translucent blue rectangle. Base centers
  inside it are selected, with live and final PyMOL synchronization.
  Shift/Ctrl adds, Alt subtracts, and Escape cancels an active rectangle.
- The ordinary cursor is an arrow. Mouse/Space panning and inertia have been
  removed. Arrows/WASD move the drawing 40 screen pixels, or 120 with Shift,
  without changing coordinates/history. Ctrl+arrows retains coordinate nudging.
  Mouse-wheel zoom, node editing, brush selection and the Gel switch remain.
- The pane's minus/close controls hide it and expand the feature sidebar;
  the top 2D button restores the same coordinates, camera and undo history.
  Hidden visual timers and transient overlays are stopped. Native window
  minimization is also enabled. Reopening returns keyboard focus to the canvas.

Pane hiding is distinct from closing the entire GUI or reanalyzing, which
still retires the current editing history. Save important layouts first.
After changing molecular coordinates or replacing a same-name object, Analyze
refreshes the GUI snapshot; CLI selection commands always analyze afresh.

## Validation

The v1.7-to-v1.8 Ubuntu PyMOL/Qt regression passes without unexpected
differences or failures. Only the explicitly requested navigation, letter
color and sequence/pane controls are normalized in the comparison. Core
feature selections, model data, layouts, brush/edit/history, persistence,
command output, PNG export and context/window cleanup remain compared.
The independent editor-batching comparison also passes, retaining v1.7's
geometry/synchronization operation counts and cached GUI/fresh CLI behavior.

New checks exercise exact rectangle membership in both gel modes, reversed
drag, add/subtract/empty selection, blue overlay cleanup and live 3D mapping;
all arrow/WASD directions and fast movement; letter formats and neutral fills;
sequence ranges, chain/residue mapping; native minimize and pane hide/restore;
camera/history preservation and keyboard focus. Synthetic residues cover T,
unknown letters, negative numbers and insertion codes.

An isolated real-Qt navigation check passes 24 Fit/scene shapes and 1,280
repeated moves at multiple zoom factors, with exact pixel steps. A sequence
component check preserves all 140 symbols including lower-case/modified bases,
chain breaks and real identifiers, and tests colors, scrolling and range input.

Qt event-interface tests do not certify physical OS mouse delivery. QtSvg is
not available in this Ubuntu PyMOL wrapper, so SVG code is retained but not
exercised. Other PyMOL builds and arbitrary structures are not certified.
PyMOL/DSSR executables are not bundled; original BSD and NAView Apache-2.0
notices remain intact. This is a local extension, not an upstream release.
