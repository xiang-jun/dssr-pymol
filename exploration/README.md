# Experimental Prototypes & Exploration (`exploration/`)

This directory houses experimental prototypes and architectural explorations that operate outside the core production release track.

---

### Directory Contents

* **`DSSR-PyMOL-GUI-v1.3.3/` (Historical Snapshot)**
  * **Lead Developer**: Eric (Enzhi) Chen
  * **Line Count**: ~10,890 lines
  * **Role**: Preserved intact as an immutable reference archive of Eric's original August 2026 distribution bundle.

* **`v2.0.0-dev/` (Active Development)**
  * **Role**: Working copy derived from the v1.3.3 bundle, subject to ongoing development, refactoring, bug fixes, and feature enhancement.
  * **Key Features**: Pure-Python RNA 2D studio (`dssr_2d`, translating features from the Jmol VARNA-plugin into PyMOL), extended GUI panels, bidirectional 2D/3D residue selection, and property inspection.

---

### Rationale & Boundaries

1. **Parallel Development Track**:
   Given the significant expansion in scope and known edge cases, developmental work in `v2.0.0-dev/` proceeds in parallel without blocking the stable `v1.1.0` production line at the repository root.
2. **Archival Integrity**:
   The `DSSR-PyMOL-GUI-v1.3.3/` subfolder remains untouched to document the exact historical baseline for auditing and technical comparison.
3. **Selective Integration**:
   Validated algorithms and isolated components from `v2.0.0-dev/` may eventually be harvested, formatted, and merged into the root plugin via standard pull requests.

---

### Usage Note

Code in this folder is intended for **testing, development, and exploration**. For production workflows and daily research use, use the verified release in the repository root.
