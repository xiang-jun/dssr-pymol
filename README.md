# dssr_select: RNA/DNA Structural Feature Selection for PyMOL

The **dssr_select** plugin is a powerful tool designed to bridge high-quality structural analysis with 3D visualization. It allows users to leverage [DSSR (Dissecting the Spatial Structure of RNA)](http://x3dna.org/) derived structural features—such as base pairs, stems, hairpins, pseudoknots, and A-minors—directly within [PyMOL](https://pymol.org/).

### ✨ Key Features

* **Interactive Graphical User Interface (GUI)**: A comprehensive Qt-based interface for browsing, filtering, and highlighting structural features without typing commands.
* **Automated Selection**: Instantly select complex motifs like junctions, multiplets, and coaxially stacked stems using DSSR’s JSON output.
* **Pseudoknot Visualization**: Identify and color individual pseudoknot layers based on dot-bracket notation.
* **Base Blocks**: Integration of the `dssr_block` logic to create schematic rectangular "block" representations of nucleic acid bases.
* **Sequence Extraction**: Quick FASTA-formatted sequence extraction with reverse-complement capabilities.

---

### 📅 Development Timeline

This plugin is the result of a multi-stage collaborative effort initiated and coordinated by **Xiang-Jun Lu**. The development followed an iterative process that merged the strengths of two researchers:

#### Phase 1: Foundational Logic (Bener Dulger)

* **Dec 06, 2025**: Initial project kick-off; development of basic nucleotide ID parsing.
* **Jan 04, 2026**: Implementation of base-pair selection and standard secondary structure detection.
* **Feb 15, 2026**: Final foundational version; added robust dot-bracket parsing for pseudoknot layers.

#### Phase 2: Feature Expansion & GUI Integration (Eric Chen)

* **Feb 18, 2026**: Expansion of the structural `FEATURE_MAP` and refined selection logic.
* **Mar 04, 2026**: Major release; introduction of the **Interactive Qt GUI**, full integration of **`dssr_block`** logic, and automated RNA structure reports.

#### Phase 3: Consolidation (Present)

* **Mar 2026**: Migration to a unified `dssr_select.py` codebase for a synergistic development environment.

---

### 📂 Repository Structure

* `dssr_select.py`: The primary, consolidated script containing all features.
* `LICENSE`: The BSD 2-Clause License.
* `NOTICE`: Detailed project history and third-party credits.
* `/archive`: Preservation of historical versions for educational and transparency purposes.
* `/bener_versions`: Original scripts from Phase 1.
* `/eric_versions`: Original scripts from Phase 2.
* `/original_utilities`: Foundational logic including the original `dssr_block.py` by Thomas Holder.

---

### 🚀 Installation & Usage

#### Requirements

1. **DSSR**: The `x3dna-dssr` executable must be installed and in your system PATH.
2. **PyMOL**: A version of PyMOL that supports `PyMOL.Qt` (standard in most Schrödinger distributions).

## Getting Started

### Installation
1. **Download**: Save the `dssr_select.py` file from this repository to your computer.
2. **Install**: In PyMOL, navigate to `Plugin` -> `Plugin Manager`.
3. **Load**: Select the `Install New Plugin` tab, click `Choose file...`, and select the `dssr_select.py` file you just downloaded.
4. **Access**: Once installed, a new **DSSR** menu item will appear under the PyMOL `Plugin` menu for easy access in future sessions.

### ⚡ Quick Demo (1-Minute Visualization)
To see the DSSR-PyMOL integration in action:

1. **Fetch a structure**: In the PyMOL console, type `fetch 1ehz`.
2. **Launch the GUI**: Go to `Plugin` -> `DSSR` (or type `dssr_gui`).
3. **Generate Blocks**: In the DSSR GUI window, click the **"make blocks"** button.

You will immediately see a schematic representation of the tRNA structure with stylized base blocks, providing a clear view of the RNA architecture.

---

### 📜 How to Cite

If you use this plugin in your research, please cite it as follows:

**Software Citation**

> Eric Chen, Bener Dulger, and Xiang-Jun Lu. **dssr_select: A PyMOL plugin for interactive RNA/DNA structural feature selection and visualization.** (2026). Available at: https://github.com/xiang-jun/dssr-pymol

**Core Technology Citation**

> Lu XJ, Bussemaker HJ, Olson WK (2015). **DSSR: an integrated software tool for dissecting the spatial structure of RNA.** *Nucleic Acids Research*, 43(21), e142.

---

### ⚖️ License

This project is licensed under the **BSD 2-Clause License**. See the [LICENSE](./LICENSE) file for the full text. This plugin incorporates components from the `dssr_block` plugin by **Thomas Holder** (c) Schrödinger LLC.
