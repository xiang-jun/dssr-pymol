# Contributing to dssr_select

Welcome! This project was built through the incredible collaborative efforts of **Bener Dulger** and **Eric Chen**. We invite the original team and the wider structural biology community to help evolve this tool.

Whether you are here to fix a bug or add a new structural feature, your contributions are what keep dssr_select synergistic and robust.

### 🌟 Our Development Philosophy
This project is built on **synergy**. We value the integration of robust structural parsing with intuitive visualization. When contributing, please consider how your changes affect both the underlying logic and the User Interface (GUI).

### 🚀 How to Contribute

#### 1. Reporting Bugs
* Check the [Issues](https://github.com/xiang-jun/dssr-pymol/issues) tab to see if the bug has already been reported.
* Use the **Bug Report** template to provide your PyMOL and DSSR versions.
* Provide a PDB ID and the specific steps to reproduce the error.

#### 2. Suggesting Enhancements
* We welcome ideas for new structural features (e.g., mapping more DSSR JSON keys to the GUI).
* Open an Issue using the **Feature Request** template to discuss the idea before starting the code.

#### 3. Submitting Code Changes
If you are a collaborator or have a fix ready:
1. **Branching**: Create a new branch for your feature (e.g., `feature/v1.1-new-motifs`).
2. **Coding Style**:
   - Maintain the existing indentation and naming conventions.
   - Ensure the `_DSSRGuiDialog` class remains decoupled from the core parsing logic where possible.
3. **Testing**: Test your changes with standard RNA structures (e.g., 1EHZ) and complex structures with pseudoknots.
4. **Pull Requests**: Submit a Pull Request (PR) to the `main` branch. Provide a clear description of what was changed and why.

### 📂 The Archive Policy
The `/archive` directory is a historical record of the project's evolution (Dec 2025 – Mar 2026).
* **Do not modify files in the `/archive` folder.** * All new development should occur in the root `dssr_select.py` file or new utility scripts in the root.

### 📜 Attribution
By contributing, you agree that your contributions will be licensed under the project's [BSD 2-Clause License](./LICENSE). New major contributors will be added to the **NOTICE** file.
