# DSSR-PyMOL GUI v1.3.3

A single-file PyMOL plugin for interactive RNA/DNA structural annotation with
X3DNA-DSSR. This research-preview build remains compatible with the upstream
`dssr-pymol` workflow and adds a pure-Python RNA 2D studio, property queries,
JSON inspection, and bidirectional 2D/3D residue selection.

This release uses a unified white scientific interface for both the main GUI
and the RNA 2D studio, including when gel rendering is enabled.

Upstream project: <https://github.com/xiang-jun/dssr-pymol>

## Package contents

- `dssr_select.py` — the complete plugin; no additional Python package is required.
- `demo_1ehz.pml` — a minimal PyMOL smoke-test session.
- `README_CN.md` — Chinese quick-start instructions.
- `RELEASE_NOTES.md` — features and validation information.
- `NOTICE` and `LICENSE` — authorship, third-party credits, and BSD-2-Clause terms.
- `SHA256SUMS.txt` — file-integrity checksums.

## Requirements

- PyMOL with Qt support (`PyMOL.Qt`).
- X3DNA-DSSR installed separately as `x3dna-dssr`.
- Network access is only needed for commands such as `fetch 1ehz`.

Jmol and Java are **not required**. The 2D viewer/editor is implemented inside
the Python plugin.

X3DNA-DSSR is separately licensed software and is intentionally not included in
this archive. Academic users can request DSSR Basic or DSSR Pro from Columbia
Technology Ventures:
<https://inventions.techventures.columbia.edu/technologies/CU20391>

## Recommended installation

1. Extract this ZIP archive.
2. Install X3DNA-DSSR and either:
   - place `x3dna-dssr` on the system `PATH`; or
   - enter its absolute path in the GUI's **exe** field.
3. Open PyMOL and choose **Plugin → Plugin Manager → Install New Plugin**.
4. Select `dssr_select.py` from this folder.
5. Restart PyMOL.

The plugin can also be loaded for one session from the PyMOL command line:

```text
run /absolute/path/to/dssr_select.py
```

## One-minute test

Enter the following commands in PyMOL:

```text
fetch 1ehz, async=0
dssr_gui
```

The DSSR window should open. Use **RNA 2D studio** in the GUI, or run:

```text
dssr_2d selection=1ehz
```

If DSSR is not on `PATH`, provide the executable explicitly:

```text
dssr_2d selection=1ehz, exe=/absolute/path/to/x3dna-dssr
```

Windows example for the GUI **exe** field:

```text
C:/Tools/DSSR/x3dna-dssr.exe
```

## Main commands

- `dssr_gui` — open the graphical interface.
- `dssr_2d selection=object_name` — open the pure-Python 2D studio.
- `dssr_select` — select a DSSR structural feature in PyMOL.
- `dssr_block` — generate DSSR base-block representations.
- `dssr_seq` — extract nucleic-acid sequence information.

## Citation

If this plugin is used in research, please cite:

Eric Chen, Bener Dulger, and Xiang-Jun Lu. *dssr_select: A PyMOL plugin for
interactive RNA/DNA structural feature selection and visualization* (2026).
<https://github.com/xiang-jun/dssr-pymol>

Please also cite the DSSR method:

Lu XJ, Bussemaker HJ, Olson WK (2015). DSSR: an integrated software tool for
dissecting the spatial structure of RNA. *Nucleic Acids Research* 43(21):e142.

## Support information

When reporting a problem, include the PyMOL version, operating system, DSSR
version, console error text, and a reproducible PDB identifier or structure.
