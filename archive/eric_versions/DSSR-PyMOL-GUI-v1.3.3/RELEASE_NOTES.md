# Release notes

## v1.3.3-white-ui-2026-08-18

- Replaces the dark-blue main GUI theme with a white scientific interface.
- Keeps the RNA 2D studio canvas pure white in both gel and flat modes.
- Retunes text, borders, disabled controls, selections, tabs, and tooltips for
  clear contrast on white backgrounds.
- Preserves gel base rendering, free brush selection, and bidirectional 2D/3D
  synchronization.

## v1.3.2-streamlined-2026-08-17

- Removes the experimental isolated-pair, kissing-loop, and H-bond feature
  categories from the command interface and GUI.
- Removes the H-bond residue, atom, and distance selection modes and their
  associated distance-object behavior.
- Keeps the remaining upstream structural-feature workflow, JSON reports and
  queries, base blocks, sequence tools, and pure-Python 2D/3D interaction.

## v1.3.1-optimized-2026-08-17

- Retains the upstream DSSR selection, sequence, base-block, and Qt GUI workflow.
- Adds DSSR property queries and report/details/full-JSON views.
- Adds save, load, copy, and lazy rendering of DSSR JSON.
- Adds a pure-Python RNA 2D studio with standard scientific layout.
- Adds optional gel-style base rendering, free brush selection, and 2D-to-3D sync.
- Adds optional 3D-to-2D synchronization and selection highlighting.
- Avoids repeated JSON formatting, coordinate-heavy selection polling, and
  unnecessary highlight-object reconstruction.
- Requires no Jmol or Java runtime.

## Validation

This release was checked with PyMOL and DSSR using PDB entry `1ehz`:

- Plugin import and command registration succeeded.
- The main `dssr_gui` window opened successfully.
- The 2D studio opened with 76 nucleotides, one chain, 21 secondary pairs, and
  13 additional DSSR pairs.
- Feature normalization, lazy JSON rendering and cache reuse, reverse 3D
  selection, and unchanged-selection highlight caching passed.
- The streamlined GUI exposes 19 structural-feature buttons and contains no
  H-bond selection-mode control.

Build date: 2026-08-18
