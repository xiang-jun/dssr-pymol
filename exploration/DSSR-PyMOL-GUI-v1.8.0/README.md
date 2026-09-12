# DSSR RNA studio 1.8.0

A single-file PyMOL plugin with one shared workspace: structural features on
the left and an editable RNA 2D drawing on the right. White background,
optional gel rendering, a linked 1D sequence, rectangle/brush selection and
bidirectional 2D/3D sync. A/C/G/T/U letters have distinct colors; node fills
and the white background stay neutral.

This is a local research extension of https://github.com/xiang-jun/dssr-pymol,
not an official upstream 1.8 release. See RELEASE_NOTES.md for this release.

## Install and run

1. Install PyMOL with `pymol.Qt` support.
2. Install X3DNA-DSSR separately. Put `x3dna-dssr` on PATH, or enter its full
   path under **Settings** in the plugin.
3. In PyMOL, use **Plugin → Plugin Manager → Install New Plugin**, select
   `dssr_select.py`, then restart PyMOL. For a temporary session, use
   `run /absolute/path/to/dssr_select.py`.
4. Load a local PDB/CIF file, or run:

```text
fetch 1ehz, async=0
dssr_gui
```

The first opening analyzes the loaded molecule. Choose a different object,
selection or state and click **Analyze** to update both panels. Changing the
analysis context retires the old drawing and its editing history. Save a
custom layout through **File → Save layout** before reanalyzing or closing.

`dssr_2d selection=1ehz` uses this same window; it does not open a second GUI.
Click an entry in the left list to select its bases in PyMOL and the 2D view.
Double-click creates a named feature selection. Use the text field to filter
entries; **Details** and **Summary** stay in the sidebar.
Selections reuse the current analysis. After editing molecular coordinates
or replacing an object under the same name, click **Analyze** to refresh it.
Console `dssr` and `dssr_select` still perform a fresh analysis on every call.

## Sequence and 2D controls

- Left-drag blank canvas to draw a blue rectangle and select base centers
  inside it, including in 3D. **Shift/Ctrl** adds; **Alt** subtracts.
- Drag a base to edit coordinates; choose a drag group from the dropdown.
- **Arrows / WASD** move the whole drawing on screen; **Shift** moves faster.
  They do not edit nucleotide coordinates. **Ctrl+arrows** nudge selected bases.
- The canvas mouse wheel zooms; **F** fits the drawing. Mouse/Space panning
  has been removed so ordinary left-drag is unambiguously a selection.
- **B** switches to freehand brush selection; **P** returns to rectangle/edit.
- The 1D strip above the canvas shows original base letters, chain breaks
  and residue IDs. Click/drag a range, **Shift** extends, **Ctrl** toggles;
  its scrollbar or wheel scrolls the sequence horizontally.
- **Ctrl+Z / Ctrl+Y** undo and redo; **Esc** clears the selection.
- **Gel** toggles both the appearance and elastic motion.
- **File** contains image export, layout save/load and DBN copying.
- **Options** contains numbering, extra pairs, letter colors and 3D sync settings.
- **− / ×** hides the sequence/2D pane, letting the feature browser fill the
  workspace. Click **2D** in the top toolbar to restore the same edited layout,
  zoom and undo history. The window itself can also be minimized.

Hiding the pane is not closing the entire GUI: closing the entire window or
reanalyzing still retires the current editing history. Save important edits.

The runtime is only `dssr_select.py`. No Java, Jmol, browser or network is needed
for local analysis. Fetching molecules requires network access. SVG export
requires QtSvg; PNG export is also supported. No DSSR executable is bundled.

## Attribution and licensing

Eric Chen, Bener Dulger, and Xiang-Jun Lu. *dssr_select: A PyMOL plugin for
interactive RNA/DNA structural feature selection and visualization* (2026).
Upstream project: https://github.com/xiang-jun/dssr-pymol

Lu XJ, Bussemaker HJ, Olson WK (2015). *DSSR: an integrated software tool for
dissecting the spatial structure of RNA*. Nucleic Acids Research 43(21):e142.

Plugin: BSD-2-Clause. The bundled NAView adaptation has separate Apache-2.0
terms, included in the source and APACHE-2.0.txt. See NOTICE. PyMOL and DSSR
are external software with their own licensing terms.
DSSR: https://inventions.techventures.columbia.edu/technologies/CU20391
