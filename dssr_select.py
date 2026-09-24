# dssr_select.py
# DSSR-based selection and 2D visualization of RNA structural features in PyMOL
#
# (c) 2026 Eric Chen, Bener Dulger, and Xiang-Jun Lu
#
# Project conceived, directed, and actively co-developed by Xiang-Jun Lu, including
#       core architecture, ongoing refactoring, bug fixes, and feature integration.
#
# CONTRIBUTIONS:
# - Eric Chen: Developed the Qt GUI ('dssr_gui'), integrated 'dssr_block', and
#              built the interactive RNA 2D layout studio ('dssr_2d').
# - Bener Dulger: Initial structural feature selection and JSON parsing; ongoing
#                 feature development, architectural refactoring, and documentation.
# - Thomas Holder: Original 'dssr_block' logic (c) Schrodinger LLC.
#
# LICENSE: BSD 2-Clause
#
# Incorporates code from 'dssr_block' and the ViennaRNA/fornac NAView engine.

from pymol import cmd, CmdException
from pymol.Qt import QtWidgets, QtCore
from pymol.plugins import addmenuitemqt
import subprocess
import json
import re
import tempfile
import os
import math
import bisect
import time
from collections import deque
from pymol.Qt import QtGui

try:
    from pymol.Qt import QtSvg
except ImportError:
    QtSvg = None

__DSSR_PLUGIN_VERSION__ = "v2.0.0-dev"
DSSR_TIMEOUT_SECONDS = 300  # Default timeout in seconds (5 minutes)

_prev_dialog = globals().get("_DSSR_GUI_DIALOG")
if _prev_dialog is not None:
    try:
        _prev_dialog.close()
        _prev_dialog.deleteLater()
    except Exception:
        pass

_HEX_COLOR_CACHE = {}

_DSSR_GUI_DIALOG = None
_DSSR_BLOCK_OBJECTS = set()
_DSSR_SELECTION_OBJECTS = set()
_DSSR_DATA_CACHE = {"key": None, "data": None}

FEATURE_MAP = {
    "pairs": "pairs",
    "stems": "stems",
    "helices": "helices",
    "stacks": "stacks",
    "nonstack": "nonStack",
    "coaxstacks": "coaxStacks",
    "atom2bases": "atom2bases",
    "aminors": "Aminors",
    "splayunits": "splayUnits",
    "hairpins": "hairpins",
    "bulges": "bulges",
    "iloops": "iloops",
    "internal": "iloops",
    "junctions": "junctions",
    "sssegments": "ssSegments",
    "ssSegments": "ssSegments",
    "multiplets": "multiplets",
    "nts": "nts",
    "pseudoknot": "pseudoknot",
    "gquadruplexes": "Gtetrads",
    "uturns": "Uturns",
}

FEATURE_ORDER = [
    "pairs",
    "stems",
    "helices",
    "stacks",
    "nonstack",
    "hairpins",
    "bulges",
    "iloops",
    "junctions",
    "sssegments",
    "multiplets",
    "coaxstacks",
    "atom2bases",
    "aminors",
    "splayunits",
    "nts",
    "pseudoknot",
    "gquadruplexes",
    "uturns",
]

FEATURE_LABELS = {
    "nonstack": "non-stack",
    "coaxstacks": "coax stacks",
    "atom2bases": "atom-base",
    "aminors": "A-minor",
    "splayunits": "splay units",
    "sssegments": "ss segments",
    "gquadruplexes": "G-tetrads",
    "uturns": "U-turns",
}

# A = Red, C = Yellow/Amber, G = Green, U/T = Cyan
PYMOL_BASE_COLORS = {
    "A": {
        "text": "#b91c1c",  # Deep crimson red (high contrast on white)
        "fill": "#fee2e2",  # Soft red pastel
        "border": "#ef4444",  # Red border
    },
    "C": {
        "text": "#b45309",  # Deep warm amber/gold (readable on white)
        "fill": "#fef9c3",  # Soft yellow pastel
        "border": "#eab308",  # Amber/yellow border
    },
    "G": {
        "text": "#15803d",  # Forest green
        "fill": "#dcfce7",  # Soft green pastel
        "border": "#22c55e",  # Green border
    },
    "U": {
        "text": "#0369a1",  # Deep cyan / sky blue
        "fill": "#e0f2fe",  # Soft cyan pastel
        "border": "#0ea5e9",  # Cyan border
    },
    "T": {
        "text": "#0369a1",  # Deep cyan / sky blue
        "fill": "#e0f2fe",  # Soft cyan pastel
        "border": "#0ea5e9",  # Cyan border
    },
    "I": {
        "text": "#6d28d9",  # Inosine / modified: purple
        "fill": "#f3e8ff",
        "border": "#a855f7",
    },
}

BLOCK_FEATURES = [
    "face",
    "edge",
    "wc",
    "g4",
    "imotif",
    "equal",
    "minor",
    "gray",
    "fill",
    "hbond",
]

LAYOUT_CHOICES = ["standard", "circular", "linear", "legacy radiate"]

LIGHT_THEME = """
QDialog, QWidget#dssrWorkspace { background: #f1f5f9; color: #0f172a; }
QMenu { background: #ffffff; color: #0f172a; border: 1px solid #cbd5e1; }
QMenu::item:selected { background: #ffe4e6; color: #e11d48; font-weight: 600; }
QLabel, QCheckBox, QGroupBox { color: #0f172a; }
QCheckBox { spacing: 5px; }
QGroupBox {
    border: 1px solid #cbd5e1; border-radius: 8px;
    margin-top: 8px; padding-top: 10px; background: #ffffff;
    font-weight: 600; color: #0f172a;
}
QGroupBox::title {
    subcontrol-origin: margin; left: 10px; padding: 0 4px;
    color: #0369a1;
}
QPushButton, QToolButton, QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {
    color: #0f172a; background: #ffffff;
    border: 1px solid #94a3b8; border-radius: 6px;
    padding: 4px 7px; min-height: 22px; font-weight: 500;
}
QPushButton:hover, QComboBox:hover, QLineEdit:focus {
    background: #f8fafc; border-color: #0284c7;
}
QPushButton:checked {
    color: #ffffff; background: #0284c7; border-color: #0369a1;
    font-weight: 700;
}
QPushButton:disabled, QComboBox:disabled {
    color: #94a3b8; background: #f1f5f9; border-color: #e2e8f0;
}
QListWidget, QPlainTextEdit {
    color: #0f172a; background: #ffffff;
    border: 1px solid #cbd5e1; border-radius: 7px;
    selection-color: #ffffff; selection-background-color: #e11d48;
}
QComboBox QAbstractItemView {
    color: #0f172a; background: #ffffff;
    selection-color: #ffffff; selection-background-color: #e11d48;
}
QTabWidget::pane { border: 1px solid #cbd5e1; border-radius: 7px; background: #ffffff; }
QTabBar::tab {
    color: #647581; background: #e2e8f0;
    border: 1px solid #cbd5e1; padding: 6px 12px;
}
QTabBar::tab:selected { color: #0369a1; background: #ffffff; font-weight: 600; }
QToolTip { color: #0f172a; background: #ffffff; border: 1px solid #0284c7; }
QPushButton:pressed { background: #e2e8f0; }
QCheckBox::indicator { width: 15px; height: 15px; }
QLabel#studioHint { color: #64748b; font-weight: 400; }
"""

DARK_THEME = """
QDialog, QWidget#dssrWorkspace { background: #0f172a; color: #f8fafc; }
QMenu { background: #1e293b; color: #f8fafc; border: 1px solid #334155; }
QMenu::item:selected { background: #e11d48; color: #ffffff; font-weight: 600; }
QLabel, QCheckBox, QGroupBox { color: #f8fafc; }
QCheckBox { spacing: 5px; }
QGroupBox {
    border: 1px solid #334155; border-radius: 8px;
    margin-top: 8px; padding-top: 10px; background: #1e293b;
    font-weight: 600; color: #38bdf8;
}
QGroupBox::title {
    subcontrol-origin: margin; left: 10px; padding: 0 4px;
    color: #38bdf8;
}
QPushButton, QToolButton, QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {
    color: #f8fafc; background: #1e293b;
    border: 1px solid #475569; border-radius: 6px;
    padding: 4px 7px; min-height: 22px; font-weight: 500;
}
QPushButton:hover, QComboBox:hover, QLineEdit:focus {
    background: #334155; border-color: #38bdf8;
}
QPushButton:checked {
    color: #ffffff; background: #0284c7; border-color: #38bdf8;
    font-weight: 700;
}
QPushButton:disabled, QComboBox:disabled {
    color: #64748b; background: #0f172a; border-color: #1e293b;
}
QListWidget, QPlainTextEdit {
    color: #f8fafc; background: #1e293b;
    border: 1px solid #334155; border-radius: 7px;
    selection-color: #ffffff; selection-background-color: #e11d48;
}
QComboBox QAbstractItemView {
    color: #f8fafc; background: #1e293b;
    selection-color: #ffffff; selection-background-color: #e11d48;
}
QTabWidget::pane { border: 1px solid #334155; border-radius: 7px; background: #1e293b; }
QTabBar::tab {
    color: #94a3b8; background: #0f172a;
    border: 1px solid #334155; padding: 6px 12px;
}
QTabBar::tab:selected { color: #38bdf8; background: #1e293b; font-weight: 600; }
QToolTip { color: #f8fafc; background: #1e293b; border: 1px solid #38bdf8; }
QPushButton:pressed { background: #334155; }
QCheckBox::indicator { width: 15px; height: 15px; }
QLabel#studioHint { color: #94a3b8; font-weight: 400; }
"""


class DssrUtils:
    @staticmethod
    def unquote(s):
        s = str(s)
        if not s:
            return s
        if s.rstrip()[-1:] not in ('"', "'"):
            return s
        return cmd.safe_eval(s)

    @staticmethod
    def _safe_tail(s, n=500):
        try:
            s = str(s)
        except Exception:
            return ""
        if len(s) <= n:
            return s
        return s[-n:]

    @staticmethod
    def _run_dssr(args, operation="DSSR", timeout=DSSR_TIMEOUT_SECONDS):
        """Run either annotation or block generation with timeout and error handling."""
        try:
            result = subprocess.run(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            raise CmdException(
                '%s timed out after %s seconds: command "%s"'
                % (operation, timeout, " ".join(args))
            )
        except OSError:
            raise CmdException('Cannot execute exe="%s"' % args[0])

        if result.returncode:
            raise CmdException(
                "%s failed (rc=%s). stderr tail: %s"
                % (
                    operation,
                    result.returncode,
                    DssrUtils._safe_tail(result.stderr),
                )
            )
        return result.stdout, result.stderr

    @staticmethod
    def run_dssr_json(pdb_path, exe, timeout=DSSR_TIMEOUT_SECONDS):
        out, err = DssrUtils._run_dssr(
            [exe, "--json", "--u-turn", "--idstr=ebi", "-i=" + pdb_path],
            operation="DSSR JSON",
            timeout=timeout,
        )
        tail = DssrUtils._safe_tail(err)
        if not out.strip():
            raise CmdException(
                "DSSR returned empty stdout (expected JSON). stderr tail: %s" % tail
            )
        try:
            return json.loads(out)
        except ValueError:
            start, end = out.find("{"), out.rfind("}")
            head = out[:120].replace("\n", " ")
            if start < 0 or end <= start:
                raise CmdException(
                    "Failed to parse DSSR JSON (no JSON object found). stdout head: %s | stderr tail: %s"
                    % (head, tail)
                )
            try:
                return json.loads(out[start : end + 1])
            except ValueError as error:
                raise CmdException(
                    "Failed to parse DSSR JSON. stdout head: %s | stderr tail: %s | err: %s"
                    % (head, tail, error)
                )

    @staticmethod
    def _selection_json(selection, state, exe, precolor=False):
        """Export and analyze a selection, always removing its temporary PDB."""
        with tempfile.TemporaryDirectory(prefix="dssr_") as directory:
            path = os.path.join(directory, "input.pdb")
            cmd.save(path, selection, state)
            if precolor:
                cmd.color("gray", selection)
            return DssrUtils.run_dssr_json(path, exe)

    @staticmethod
    def _invalidate_cache():
        """Reset the shared DSSR analysis cache."""
        global _DSSR_DATA_CACHE
        _DSSR_DATA_CACHE["key"] = None
        _DSSR_DATA_CACHE["data"] = None

    @staticmethod
    def _cached_selection_json(selection, state, exe, precolor=False, force=False):
        """Analyze a selection with DSSR, reusing parsed JSON if context and atom count match."""
        global _DSSR_DATA_CACHE

        try:
            current_count = int(cmd.count_atoms(selection, state=state))
        except Exception:
            current_count = -1

        cache_key = (str(selection), int(state), str(exe), current_count)

        if (
            not force
            and _DSSR_DATA_CACHE["key"] == cache_key
            and _DSSR_DATA_CACHE["data"] is not None
        ):
            if precolor:
                cmd.color("gray", selection)
            return _DSSR_DATA_CACHE["data"]

        # Run fresh analysis and update the shared cache
        data = DssrUtils._selection_json(selection, state, exe, precolor=precolor)
        _DSSR_DATA_CACHE["key"] = cache_key
        _DSSR_DATA_CACHE["data"] = data
        return data

    @staticmethod
    def _hex_to_rgb01(h):

        s = str(h).strip()

        if s.startswith('"') and s.endswith('"'):
            s = s[1:-1].strip()
        if s.startswith("'") and s.endswith("'"):
            s = s[1:-1].strip()

        if s.lower().startswith("0x"):
            s = s[2:]
        if s.startswith("#"):
            s = s[1:]

        s = s.strip()
        if re.fullmatch(r"[0-9a-fA-F]{3}", s):
            s = "".join([c * 2 for c in s])

        if not re.fullmatch(r"[0-9a-fA-F]{6}", s):
            raise CmdException(
                'Invalid hex color "%s". Use FF00AA or 0xFF00AA (or "#FF00AA" quoted).'
                % h
            )

        r = int(s[0:2], 16) / 255.0
        g = int(s[2:4], 16) / 255.0
        b = int(s[4:6], 16) / 255.0
        return s.lower(), [r, g, b]

    @staticmethod
    def _resolve_color_spec(color_spec):
        if color_spec is None:
            return None
        s = str(color_spec).strip()
        if not s:
            return None
        if s.lower() in ("auto", "default"):
            return None

        is_hexish = (
            s.startswith("#")
            or s.lower().startswith("0x")
            or (len(s) in (3, 6) and all(c in "0123456789abcdefABCDEF" for c in s))
        )
        if is_hexish:
            hex6, rgb = DssrUtils._hex_to_rgb01(s)
            if hex6 in _HEX_COLOR_CACHE:
                return _HEX_COLOR_CACHE[hex6]
            cname = "dssr_hex_%s" % hex6
            cmd.set_color(cname, rgb)
            _HEX_COLOR_CACHE[hex6] = cname
            return cname

        name = s.lower()
        try:
            idx = cmd.get_color_index(name)
            if idx < 0:
                raise Exception("unknown")
        except Exception:
            raise CmdException(
                'Unknown color "%s". Use a PyMOL color name (e.g., blue) or hex like FF00AA / 0xFF00AA.'
                % s
            )
        return name

    @staticmethod
    def matches_boolean_query(text, query):
        """Evaluate whether `text` satisfies a boolean query supporting:
        - AND: space or 'and' / '&&'
        - OR: 'or' / '|'
        - NOT: '-' / '!' / 'not'
        """
        text = str(text).lower()
        query = str(query).strip().lower()
        if not query:
            return True

        # 1. Split across OR clauses (either '|' or whole-word 'or')
        or_clauses = [
            c.strip() for c in re.split(r"\s+\bor\b\s+|\|", query) if c.strip()
        ]

        # 2. Entry matches if it satisfies ANY OR clause
        for clause in or_clauses:
            tokens = [t.strip() for t in clause.split() if t.strip()]
            clause_match = True
            expect_not = False

            for token in tokens:
                if token in ("and", "&&"):
                    continue
                if token in ("not", "!"):
                    expect_not = True
                    continue

                # Handle negative terms like "-wc", "!wobble", or "NOT wc"
                if expect_not or token.startswith(("-", "!")):
                    neg_term = token.lstrip("-!") if not expect_not else token
                    expect_not = False
                    if neg_term and neg_term in text:
                        clause_match = False
                        break
                else:
                    # Positive term must appear in the entry text
                    if token not in text:
                        clause_match = False
                        break

            if clause_match:
                return True

        return False

    @staticmethod
    def base_style(base, enabled=True):
        b = str(base).strip().upper()
        if not enabled:
            return {
                "text": QtGui.QColor("#0f172a"),
                "fill": QtGui.QColor("#ffffff"),
                "border": QtGui.QColor(30, 41, 59),
            }
        if b in PYMOL_BASE_COLORS:
            spec = PYMOL_BASE_COLORS[b]
        elif b in ("P", "PSU"):  # Pseudouridine
            spec = PYMOL_BASE_COLORS["U"]
        else:
            spec = PYMOL_BASE_COLORS.get(
                "I",
                {
                    "text": "#6d28d9",
                    "fill": "#f3e8ff",
                    "border": "#a855f7",
                },
            )
        return {
            "text": QtGui.QColor(spec["text"]),
            "fill": QtGui.QColor(spec["fill"]),
            "border": QtGui.QColor(spec["border"]),
        }

    @staticmethod
    def base_text_color(base, enabled=True):
        return DssrUtils.base_style(base, enabled)["text"]


class DssrParser:
    @staticmethod
    def feature_entries(dssr_data, feature):
        """Return a normalized list for a FEATURE_MAP entry.

        DSSR normally emits arrays, but ``nonStack`` can be emitted as a
        single object. Keeping that compatibility rule here prevents the GUI,
        command API, query engine, and selection builder from drifting apart.
        """
        if not isinstance(dssr_data, dict):
            return []
        feature = str(feature).strip()
        json_key = FEATURE_MAP.get(feature)
        if not json_key:
            return []
        value = dssr_data.get(json_key)
        if feature == "nonstack" and isinstance(value, dict):
            if value.get("num_nts", 0) > 0 or value.get("nts_long"):
                return [value]
            return []
        return value if isinstance(value, list) else []

    @staticmethod
    def parse_nt_id(nt_id):
        """
        Parses a nucleotide identifier in Jmol/EBI Unit ID format:
        |Model Number| Chain ID|Residue Identifier|Residue Number|Atom Name|Alternate ID|Insertion Code|
        Returns a compatible tuple of (chain, resi_number).
        """
        parts = str(nt_id).strip().split("|")
        if len(parts) >= 5:
            return parts[2].strip(), parts[4].strip()
        raise CmdException('Unexpected nt_id format: "%s"' % nt_id)

    @staticmethod
    def parse_a2b_atom(atom_id):
        """
        Parses atom-to-base strings in Jmol/EBI Unit ID format.
        Returns a compatible tuple of (chain, resi_number, atom_name).
        """
        parts = str(atom_id).strip().split("|")
        if len(parts) >= 6:
            return parts[2].strip(), parts[4].strip(), parts[5].strip()
        raise CmdException('Unexpected atom format: "%s"' % atom_id)

    @staticmethod
    def parse_dotbracket_pseudoknots(dotbracket):
        """Return noncanonical bracket/letter layers; separators consume no index."""
        openers = "([{<"
        closer_to_open = dict(
            zip(
                ")]}>abcdefghijklmnopqrstuvwxyz", openers + "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            )
        )
        layer_for = dict(zip(openers, range(4)))
        stacks, layers = {}, {}
        nt_index = 0
        for char in str(dotbracket):
            if char.isspace() or char in "&+":
                continue
            if char in openers or "A" <= char <= "Z":
                stacks.setdefault(char, []).append(nt_index)
                if char not in layer_for:
                    layer_for[char] = len(layer_for)
            elif char in closer_to_open:
                opener = closer_to_open[char]
                stack = stacks.get(opener)
                if stack:
                    start = stack.pop()
                    layer = layer_for[opener]
                    if layer:
                        layers.setdefault(layer, []).append((start, nt_index))
            nt_index += 1
        for pairs in layers.values():
            pairs.sort()
        return layers

    @staticmethod
    def _sort_resi_key(resi_str):
        """Sort residue numbers naturally (numeric when possible, fallback to string)."""
        try:
            return (0, int(str(resi_str)))
        except Exception:
            return (1, str(resi_str))

    @staticmethod
    def _compact_sel_from_residues(residues):
        """Build a compact PyMOL selection string grouping residues by chain with '+'.

        Example:
            Input:  {('A', '10'), ('A', '11'), ('B', '5')}
            Output: '(chain A and resi 10+11) or (chain B and resi 5)'
        """
        if not residues:
            return ""

        by_chain = {}
        no_chain = set()

        for item in residues:
            if isinstance(item, (tuple, list)) and len(item) >= 2:
                c = str(item[0]).strip()
                r = str(item[1]).strip()
            else:
                continue

            if not r:
                continue

            if c:
                by_chain.setdefault(c, set()).add(r)
            else:
                no_chain.add(r)

        parts = []
        for c in sorted(by_chain.keys()):
            resis = sorted(by_chain[c], key=DssrParser._sort_resi_key)
            resi_expr = "+".join(resis)
            parts.append("(chain %s and resi %s)" % (c, resi_expr))

        if no_chain:
            resis = sorted(no_chain, key=DssrParser._sort_resi_key)
            parts.append("(resi %s)" % "+".join(resis))

        return " or ".join(parts)

    @staticmethod
    def _residue_selection(residues):
        """Direct all residue selection generation through the compact generator."""
        return DssrParser._compact_sel_from_residues(residues)

    @staticmethod
    def build_selection_from_layer(layer_pairs, nts_list):
        residues = set()
        for open_index, close_index in layer_pairs:
            for index in (open_index, close_index):
                if index < len(nts_list) and nts_list[index].get("nt_id"):
                    residues.add(DssrParser.parse_nt_id(nts_list[index]["nt_id"]))
        return DssrParser._residue_selection(residues) or None

    @staticmethod
    def build_selection_from_pair(pair_entry):
        nt1 = pair_entry.get("nt1")
        nt2 = pair_entry.get("nt2")
        if not nt1 or not nt2:
            raise CmdException("Pair entry missing nt1 or nt2")
        residues = {DssrParser.parse_nt_id(nt1), DssrParser.parse_nt_id(nt2)}
        return DssrParser._residue_selection(residues)

    @staticmethod
    def build_selection_from_nts_list(nts_list):
        if not nts_list:
            raise CmdException("Empty nucleotide list")
        residues = {DssrParser.parse_nt_id(nt) for nt in nts_list}
        return DssrParser._residue_selection(residues)

    @staticmethod
    def _pair_residues(pairs):
        return {
            DssrParser.parse_nt_id(pair[key])
            for pair in pairs
            for key in ("nt1", "nt2")
            if pair.get(key)
        }

    @staticmethod
    def build_selection_from_stem(stem_entry):
        pairs = stem_entry.get("pairs", [])
        if not pairs:
            raise CmdException("Stem has no pairs")
        return DssrParser._residue_selection(DssrParser._pair_residues(pairs))

    @staticmethod
    def build_selection_from_hairpin(hairpin_entry):
        nts_long = hairpin_entry.get("nts_long")
        if not nts_long:
            raise CmdException("Hairpin missing nts_long field")
        nts_list = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
        return DssrParser.build_selection_from_nts_list(nts_list)

    @staticmethod
    def build_selection_from_coaxstack(coax_entry, stems_list):
        stem_indices = coax_entry.get("stem_indices", [])
        if not stem_indices:
            raise CmdException("coaxStacks entry missing stem_indices")
        residues = set()
        for index in stem_indices:
            try:
                index = int(index)
            except (TypeError, ValueError, OverflowError):
                continue
            if 1 <= index <= len(stems_list):
                residues.update(
                    DssrParser._pair_residues(stems_list[index - 1].get("pairs", []))
                )
        if not residues:
            raise CmdException("Could not build selection for coaxStacks entry")
        return DssrParser._residue_selection(residues)

    @staticmethod
    def build_selection_from_atom2base(a2b_entry):
        atom = a2b_entry.get("atom")
        nt = a2b_entry.get("nt")

        clauses = []

        if nt:
            c_nt, r_nt = DssrParser.parse_nt_id(nt)
            clauses.append(DssrParser._compact_sel_from_residues([(c_nt, r_nt)]))

        if atom:
            c_a, r_a, atom_name = DssrParser.parse_a2b_atom(atom)
            atom_name = str(atom_name).replace('"', '\\"')
            if c_a:
                clauses.append(
                    '(chain %s and resi %s and name "%s")' % (c_a, r_a, atom_name)
                )
            else:
                clauses.append('(resi %s and name "%s")' % (r_a, atom_name))

        if not clauses:
            raise CmdException("atom2bases entry missing atom and nt")
        return " or ".join(clauses)

    @staticmethod
    def build_selection_from_aminor(aminor_entry):
        description = aminor_entry.get("desc_long", "")
        if not description or "vs" not in description:
            raise CmdException("Aminors entry missing desc_long")
        left, right = description.split("vs", 1)
        nts = [nt.strip() for nt in [left] + right.split(",") if nt.strip()]
        if not nts:
            raise CmdException("Aminors entry has empty residues")
        return DssrParser.build_selection_from_nts_list(nts)

    @staticmethod
    def build_selection_from_gquadruplex(gquad_entry):
        nts_long = gquad_entry.get("nts_long", "")
        nts_list = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
        return DssrParser.build_selection_from_nts_list(nts_list)

    # Both annotations store their residues in the same nts_long field.
    build_selection_from_uturns = build_selection_from_gquadruplex

    @staticmethod
    def _shorten_nts_long(nts_long, max_items=6):
        if not nts_long:
            return ""
        nts = [x.strip() for x in nts_long.split(",") if x.strip()]
        if len(nts) <= max_items:
            return ", ".join(nts)
        half = max_items // 2
        return "%s, ..., %s" % (", ".join(nts[:half]), ", ".join(nts[-half:]))

    @staticmethod
    def _preview_entry(feature, entry, i):
        if feature == "pairs":
            nt1 = entry.get("nt1", "?")
            nt2 = entry.get("nt2", "?")
            lw = entry.get("LW", entry.get("bp", ""))
            name = entry.get("name", "").strip()

            # Combine LW classification and pair name (e.g., "cWW, WC")
            tags = [t for t in (lw, name) if t]
            tag_str = (" (%s)" % ", ".join(tags)) if tags else ""
            return "%d: %s - %s%s" % (i, nt1, nt2, tag_str)

        if feature in ("stems", "helices"):
            n = (
                len(entry.get("pairs", []))
                if isinstance(entry.get("pairs", []), list)
                else 0
            )
            nm = entry.get("name", entry.get("index", ""))
            return "%d: %s (pairs=%d)" % (i, str(nm), n)

        if feature in (
            "stacks",
            "nonstack",
            "hairpins",
            "bulges",
            "iloops",
            "internal",
            "junctions",
            "sssegments",
            "ssSegments",
            "multiplets",
            "splayunits",
            "gquadruplexes",
            "uturns",
        ):
            s = DssrParser._shorten_nts_long(entry.get("nts_long", ""))
            return "%d: %s" % (i, s if s else "(missing nts_long)")

        if feature == "coaxstacks":
            return "%d: helix=%s stems=%s" % (
                i,
                str(entry.get("helix_index", "")),
                str(entry.get("stem_indices", [])),
            )

        if feature == "atom2bases":
            t = entry.get("type", "")
            atom = entry.get("atom", "?")
            nt = entry.get("nt", "?")
            return "%d: %s atom=%s nt=%s" % (i, t if t else "entry", atom, nt)

        if feature == "aminors":
            ds = entry.get("desc_short", "")
            dl = entry.get("desc_long", "")
            return "%d: %s" % (i, ds if ds else dl)

        if feature == "nts":
            return "%d: %s" % (i, entry.get("nt_id", "?"))

        return "%d: (no preview)" % i

    @staticmethod
    def _extract_dotbracket(dssr_data):
        if "dbn" not in dssr_data:
            raise CmdException("No dot-bracket notation found in DSSR output")
        dbn_data = dssr_data["dbn"]

        if isinstance(dbn_data, dict):
            if (
                isinstance(dbn_data.get("all_chains"), dict)
                and "sstr" in dbn_data["all_chains"]
            ):
                return dbn_data["all_chains"]["sstr"]
            if "sstr" in dbn_data:
                return dbn_data["sstr"]
            for v in dbn_data.values():
                if isinstance(v, dict) and "sstr" in v:
                    return v["sstr"]
            raise CmdException("Could not find sstr field in dbn data")

        return dbn_data

    @staticmethod
    def _extract_chain_names(dssr_data):
        chains = set()
        nts = dssr_data.get("nts", []) if isinstance(dssr_data, dict) else []
        if isinstance(nts, list):
            for nt in nts:
                try:
                    chain, _ = DssrParser.parse_nt_id(nt.get("nt_id", ""))
                    if chain:
                        chains.add(chain)
                except (AttributeError, CmdException):
                    continue
        return sorted(chains)

    @staticmethod
    def _count_pseudoknot_layers(dssr_data):
        try:
            dotbracket = DssrParser._extract_dotbracket(dssr_data)
            layers = DssrParser.parse_dotbracket_pseudoknots(dotbracket)
            return len(layers) if layers else 0
        except Exception:
            return 0

    @staticmethod
    def _format_rna_summary_text(dssr_data):
        chains = DssrParser._extract_chain_names(dssr_data)
        lines = [
            "RNA Structure Summary",
            "---------------------",
            "Chains: %s" % (" ".join(chains) if chains else "(unknown)"),
        ]
        for label, key in (
            ("Base pairs", "pairs"),
            ("Hairpins", "hairpins"),
            ("Stems", "stems"),
            ("Bulges", "bulges"),
            ("Junctions", "junctions"),
            ("Pseudoknots", None),
            ("A-minor interactions", "Aminors"),
            ("Stacking interactions", "stacks"),
            ("G-quadruplexes", "Gtetrads"),
            ("U-turns", "Uturns"),
        ):
            entries = dssr_data.get(key) if key else None
            count = (
                (len(entries) if isinstance(entries, list) else 0)
                if key
                else (DssrParser._count_pseudoknot_layers(dssr_data))
            )
            lines.append("%s: %d" % (label, count))
        return "\n".join(lines)

    @staticmethod
    def _build_residue_sel_from_dssr(dssr_data, feature, index):
        feature = str(feature).lower().strip()
        idx = int(index)

        if feature == "pseudoknot":
            dotbracket = DssrParser._extract_dotbracket(dssr_data)
            nts_list = dssr_data.get("nts", None)
            if nts_list is None:
                raise CmdException("No nts found in DSSR output")

            layers = DssrParser.parse_dotbracket_pseudoknots(dotbracket)
            if not layers:
                raise CmdException("No pseudoknot layers found")

            layer_keys = sorted(layers.keys())
            if idx < 1 or idx > len(layer_keys):
                raise CmdException(
                    "pseudoknot layer index %d out of range (1..%d)"
                    % (idx, len(layer_keys))
                )

            pairs = layers[layer_keys[idx - 1]]
            sel_str = DssrParser.build_selection_from_layer(pairs, nts_list)
            if not sel_str:
                raise CmdException(
                    "Could not build selection for pseudoknot layer %d" % idx
                )
            return sel_str

        if feature not in FEATURE_MAP:
            raise CmdException('Unknown feature "%s"' % feature)

        json_key = FEATURE_MAP[feature]
        feature_list = DssrParser.feature_entries(dssr_data, feature)
        if not feature_list:
            raise CmdException('No "%s" found in DSSR output' % json_key)

        if idx < 1 or idx > len(feature_list):
            raise CmdException(
                "%s index %d out of range (1..%d)" % (feature, idx, len(feature_list))
            )

        entry = feature_list[idx - 1]

        builders = {
            "pairs": DssrParser.build_selection_from_pair,
            "stems": DssrParser.build_selection_from_stem,
            "helices": DssrParser.build_selection_from_stem,
            "hairpins": DssrParser.build_selection_from_hairpin,
            "atom2bases": DssrParser.build_selection_from_atom2base,
            "aminors": DssrParser.build_selection_from_aminor,
            "gquadruplexes": DssrParser.build_selection_from_gquadruplex,
            "uturns": DssrParser.build_selection_from_uturns,
        }
        if feature in builders:
            return builders[feature](entry)

        if feature in (
            "stacks",
            "nonstack",
            "bulges",
            "iloops",
            "internal",
            "junctions",
            "sssegments",
            "ssSegments",
            "multiplets",
            "splayunits",
            "uturns",
        ):
            nts_long = entry.get("nts_long", "")
            if not nts_long:
                raise CmdException("%s entry missing nts_long field" % feature)
            nts_list_parsed = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
            return DssrParser.build_selection_from_nts_list(nts_list_parsed)

        if feature == "coaxstacks":
            stems_list = dssr_data.get("stems", [])
            if not stems_list:
                raise CmdException("No stems found, required for coaxStacks")
            return DssrParser.build_selection_from_coaxstack(entry, stems_list)

        if feature == "nts":
            nt_id = entry.get("nt_id")
            if not nt_id:
                raise CmdException("Nucleotide entry missing nt_id field")
            c, r = DssrParser.parse_nt_id(nt_id)
            return "(chain %s and resi %s)" % (c, r)

        raise CmdException('Feature "%s" not supported for residue selection' % feature)


class DssrCmd:
    @staticmethod
    def dssr_select(
        selection="all",
        state=-1,
        feature="pairs",
        index=1,
        name="dssr_select",
        exe="x3dna-dssr",
        show_info=0,
        quiet=1,
        color="auto",
        precolor=1,
    ):

        state = int(state)
        index = int(index)
        show_info = int(show_info)
        quiet = int(quiet)
        precolor = int(precolor)
        feature = DssrUtils.unquote(feature).lower().strip()

        user_color = DssrUtils._resolve_color_spec(DssrUtils.unquote(color).strip())

        if feature in ("features", "help"):
            keys = sorted(FEATURE_MAP.keys())
            print("Supported features: " + ", ".join(keys))
            print("Tip: use index=0 to list detected items for a feature.")
            return

        if feature not in FEATURE_MAP:
            valid = ", ".join(sorted(FEATURE_MAP.keys()))
            raise CmdException('Unknown feature "%s". Valid: %s' % (feature, valid))

        if state == 0 or state < 0:
            state = cmd.get_state()

        dssr_data = DssrUtils._cached_selection_json(
            selection, state, exe, precolor=bool(precolor)
        )
        return DssrCmd._select_feature_data(
            dssr_data,
            selection,
            state,
            feature,
            index,
            name,
            user_color,
            show_info,
            quiet,
        )

    @staticmethod
    def _create_feature_selection(name, selection, residue_selection, quiet=0):
        cmd.select(name, "((%s) and (%s))" % (selection, residue_selection))
        if int(cmd.count_atoms(name)) <= 0:
            cmd.delete(name)
            raise CmdException(
                "DSSR residues did not map back to the requested PyMOL selection"
            )
        _DSSR_SELECTION_OBJECTS.add(str(name))

        try:
            cmd.delete("indicate")
        except Exception:
            pass

        # Force a clean 0 -> 1 transition so PyMOL activates the pink selection indicators
        cmd.disable(name)
        cmd.enable(name)

        if not quiet:
            _sel_residues = set()
            cmd.iterate(
                name,
                "_sel_residues.add((chain, resi))",
                space={"_sel_residues": _sel_residues},
            )
            print(
                "dssr_select: %s" % DssrParser._compact_sel_from_residues(_sel_residues)
            )

    @staticmethod
    def _dssr_default_selection():
        objs = cmd.get_object_list("enabled")
        if len(objs) == 1:
            return objs[0]
        return "all"

    @staticmethod
    def _dssr(
        sel=None,
        selection=None,
        f=None,
        feature="pairs",
        i=None,
        index=1,
        n=None,
        name=None,
        q=None,
        quiet=0,
        si=None,
        show_info=0,
        st=None,
        state=-1,
        exe="x3dna-dssr",
        color="auto",
        display=0,
        stick_radius=0.25,
        do_zoom=1,
        pc=None,
        precolor=1,
    ):
        selection = selection or sel or DssrCmd._dssr_default_selection()

        feature_in = f if f is not None else feature
        feature_in = DssrUtils.unquote(feature_in).strip()

        if feature_in.lower() in ("features", "help"):
            DssrCmd.dssr_select(
                selection=selection,
                state=state,
                feature="features",
                index=0,
                name="dssr_select",
                exe=exe,
                show_info=0,
                quiet=int(q if q is not None else quiet),
                color="auto",
                precolor=int(precolor),
            )
            return

        idx = int(i if i is not None else index)
        qt = int(q if q is not None else quiet)
        si2 = int(si if si is not None else show_info)
        st2 = int(st if st is not None else state)

        if pc is not None:
            precolor = int(pc)
        precolor = int(precolor)

        nm = name if name is not None else n
        if not nm:
            nm = "%s%d" % (feature_in.lower(), idx)

        DssrCmd.dssr_select(
            selection=selection,
            state=st2,
            feature=feature_in,
            index=idx,
            name=nm,
            exe=exe,
            show_info=si2,
            quiet=qt,
            color=color,
            precolor=precolor,
        )
        DssrCmd._display_feature_selection(nm, display, stick_radius, do_zoom)

    @staticmethod
    def _unused_name(prefix):
        try:
            return cmd.get_unused_name(prefix)
        except Exception:
            base = str(prefix) if prefix else "obj"
            name = base
            k = 1
            while True:
                try:
                    exists = name in cmd.get_object_list()
                except Exception:
                    exists = False
                if not exists:
                    return name
                k += 1
                name = "%s%d" % (base, k)

    @staticmethod
    def dssr_block(
        selection="all",
        state=-1,
        block_file="face",
        block_depth=0.5,
        block_color="",
        name="",
        exe="x3dna-dssr",
        quiet=1,
    ):
        """
        DESCRIPTION

            Create a nucleic acid base "block" cartoon with DSSR.

            Requires the "x3dna-dssr" program, available from URL:
                https://inventions.techventures.columbia.edu/technologies/dssr-an-integrated-software--CU20391

        USAGE

            dssr_block [ selection [, state [, block_file [, block_depth
                [, block_color [, name [, exe ]]]]]]]

        ARGUMENTS

            selection = str: atom selection {default: all}

            state = int: object state (0 for all states) {default: -1, current state}

            block_file = face|edge|wc|g4|imotif|equal|minor|gray|fill|hbond:
                         Corresponds to the --block-file option (see DSSR manual).
                         Values can be combined, e.g. "wc-minor". {default: face}

            block_depth = float: thickness of rectangular blocks {default: 0.5}

            block_color = str: Corresponds to the --block-color option {default: }

            name = str: name of new CGO object {default: dssr_block##}

            exe = str: path to "x3dna-dssr" executable {default: x3dna-dssr}

        EXAMPLE

            fetch 1ehz, async=0
            as cartoon
            dssr_block
            set cartoon_ladder_radius, 0.1
            set cartoon_ladder_color, gray
            set cartoon_nucleic_acid_mode, 1

            # multi-state
            fetch 2n2d, async=0
            dssr_block 2n2d, 0
            set all_states

            # custom coloring
            fetch 1msy, async=0
            dssr_block block_color=N red | minor 0.9 | major yellow
        """
        try:
            state = int(state)
        except Exception:
            state = -1
        quiet = int(quiet)

        if state < 0:
            try:
                state = int(cmd.get_state())
            except Exception:
                state = 1

        if not name:
            name = DssrCmd._unused_name("dssr_block")

        with tempfile.TemporaryDirectory(prefix="dssr_block_") as directory:
            tmpfilepdb = os.path.join(directory, "input.pdb")
            tmpfiler3d = os.path.join(directory, "blocks.r3d")
            if state == 0:
                try:
                    n_states = int(cmd.count_states(selection))
                except Exception:
                    n_states = 1
                states = list(range(1, max(1, n_states) + 1))
            else:
                states = [max(1, int(state))]

            for st in states:
                cmd.save(tmpfilepdb, selection, st)

                args = [
                    exe,
                    "--block-file=" + DssrUtils.unquote(block_file),
                    "--block-depth=" + str(block_depth),
                    "-i=" + tmpfilepdb,
                    "-o=" + tmpfiler3d,
                ]

                # Incorporate block_color argument
                if block_color:
                    args.append("--block-color=" + DssrUtils.unquote(block_color))

                DssrUtils._run_dssr(
                    args, operation="DSSR block", timeout=DSSR_TIMEOUT_SECONDS
                )

                cmd.load(tmpfiler3d, name, max(1, st), zoom=0)

            _DSSR_BLOCK_OBJECTS.add(str(name))
            if not quiet:
                print(
                    'dssr_block: loaded "%s" (block_file=%s, block_depth=%s)'
                    % (name, str(block_file), str(block_depth))
                )

    @staticmethod
    def _select_feature_data(
        dssr_data,
        selection,
        state,
        feature,
        index,
        name,
        user_color=None,
        show_info=0,
        quiet=1,
    ):
        """Create a feature selection from already analyzed DSSR data."""
        json_key = FEATURE_MAP[feature]
        if feature == "pseudoknot":
            layer_colors = ["blue", "pink", "green", "yellow", "orange"]

            dotbracket = DssrParser._extract_dotbracket(dssr_data)
            nts_list = dssr_data.get("nts", None)
            if nts_list is None:
                raise CmdException("No nts found in DSSR output")

            layers = DssrParser.parse_dotbracket_pseudoknots(dotbracket)
            if not layers:
                raise CmdException("No pseudoknot layers found in structure")

            layer_keys = sorted(layers.keys())

            if index == 0:
                print("pseudoknot: %d layer(s)" % len(layer_keys))
                for j, k in enumerate(layer_keys, 1):
                    print(
                        "  layer %d (key=%s): %d pair(s)" % (j, str(k), len(layers[k]))
                    )
                if not quiet and show_info:
                    print("pseudoknot dot-bracket: " + str(dotbracket))
                return

            if index < 1 or index > len(layer_keys):
                raise CmdException(
                    "Layer index %d out of range (1..%d)" % (index, len(layer_keys))
                )

            layer_key = layer_keys[index - 1]
            pairs = layers[layer_key]
            layer_color = layer_colors[(index - 1) % len(layer_colors)]

            sel_str = DssrParser.build_selection_from_layer(pairs, nts_list)
            if sel_str is None:
                raise CmdException("Could not build selection for layer %d" % index)

            DssrCmd._create_feature_selection(name, selection, sel_str, quiet=quiet)
            cmd.color(user_color if user_color else layer_color, name)

            if not quiet:
                print(
                    'dssr_select: selection "%s" pseudoknot layer %d with %d pair(s)'
                    % (name, index, len(pairs))
                )
            return

        feature_list = DssrParser.feature_entries(dssr_data, feature)
        if not feature_list:
            raise CmdException('No "%s" found in DSSR output' % json_key)

        if index == 0:
            total = len(feature_list)
            print("%s: %d item(s)" % (feature, total))
            show_n = 20 if not quiet else 10
            show_n = min(show_n, total)
            for i in range(show_n):
                print("  " + DssrParser._preview_entry(feature, feature_list[i], i + 1))
            if total > show_n:
                print("  ... (%d more)" % (total - show_n))
            return

        if index < 1 or index > len(feature_list):
            raise CmdException(
                "%s index %d out of range (1..%d)" % (feature, index, len(feature_list))
            )

        sel_str = DssrParser._build_residue_sel_from_dssr(dssr_data, feature, index)
        if not sel_str:
            raise CmdException(
                "Could not build selection for %s index %d" % (feature, index)
            )
        DssrCmd._create_feature_selection(name, selection, sel_str, quiet=quiet)
        cmd.color(user_color if user_color else "pink", name)

        if not quiet:
            print(
                'dssr_select: created selection "%s" for %s (index %d) in state %d'
                % (name, feature, index, state)
            )

    @staticmethod
    def _display_feature_selection(name, display=0, stick_radius=0.25, do_zoom=1):
        if int(display):
            cmd.show("sticks", name)
            try:
                cmd.set("stick_radius", float(stick_radius), name)
            except Exception:
                pass
            if int(do_zoom):
                cmd.zoom(name)

    @staticmethod
    def dssr_2d(
        selection="all",
        state=-1,
        exe="x3dna-dssr",
        layout="standard",
        number_every=10,
        show_tertiary=0,
        title="",
        quiet=1,
    ):
        global _DSSR_GUI_DIALOG
        selection = DssrUtils.unquote(selection)
        exe = DssrUtils.unquote(exe)
        layout = DssrUtils.unquote(layout)
        title = DssrUtils.unquote(title)
        state, number_every = int(state), int(number_every)
        if state <= 0:
            state = int(cmd.get_state())
        if _DSSR_GUI_DIALOG is None:
            _DSSR_GUI_DIALOG = DssrGuiDialog()
        host = _DSSR_GUI_DIALOG
        if host._analysis_context != (selection, state, exe):
            host._clear_analysis("Analyzing the requested structure...")
        try:
            data = host._get_dssr_data(selection, state, exe, 0)
            editor = host.show_analysis(
                data,
                selection,
                state,
                exe,
                algorithm=layout,
                number_every=number_every,
                show_tertiary=int(show_tertiary),
                title=title,
            )
        except Exception as error:
            host._clear_analysis("Analysis error: %s" % error)
            raise
        host.show()
        host.raise_()
        host.activateWindow()
        host.show_2d_btn.setChecked(True)
        editor.view.setFocus(QtCore.Qt.OtherFocusReason)
        if not int(quiet):
            print("dssr_2d: shared workspace — %s" % editor.model.summary())
        return editor


class DssrUI:
    """Convenience factory and setup routines for Qt widgets and graphics items."""

    @staticmethod
    def button(text, clicked, tip=""):
        widget = QtWidgets.QPushButton(text)
        widget.clicked.connect(clicked)
        widget.setToolTip(tip)
        return widget

    @staticmethod
    def checkbox(text, checked=False, changed=None, tip=""):
        widget = QtWidgets.QCheckBox(text)
        widget.setChecked(checked)
        if changed is not None:
            widget.toggled.connect(changed)
        widget.setToolTip(tip)
        return widget

    @staticmethod
    def combo(items, editable=False, tip=""):
        widget = QtWidgets.QComboBox()
        widget.setEditable(editable)
        widget.addItems(items)
        widget.setToolTip(tip)
        return widget

    @staticmethod
    def spinbox(minimum, maximum, value, decimals=False, step=1, suffix="", tip=""):
        widget = QtWidgets.QDoubleSpinBox() if decimals else QtWidgets.QSpinBox()
        widget.setRange(minimum, maximum)
        widget.setSingleStep(step)
        widget.setValue(value)
        widget.setSuffix(suffix)
        widget.setToolTip(tip)
        return widget

    @staticmethod
    def menu_button(text, actions, parent):
        button = QtWidgets.QToolButton(parent)
        button.setText(text)
        button.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        menu = QtWidgets.QMenu(button)
        for label, callback in actions:
            menu.addAction(label, callback)
        button.setMenu(menu)
        return button

    @staticmethod
    def no_mouse(item):
        """Prevent item from intercepting mouse clicks, letting underlying views or parent items handle them."""
        try:
            item.setAcceptedMouseButtons(QtCore.Qt.NoButton)
        except Exception:
            pass


class DssrGuiDialog(QtWidgets.QDialog if QtWidgets else object):
    """One analysis context, feature browser, and embedded RNA editor."""

    PAGE_SIZE = 500

    def __init__(self):
        super().__init__()
        self.setWindowTitle("DSSR RNA studio")
        self.setWindowFlag(QtCore.Qt.WindowMinimizeButtonHint, True)
        self.resize(1360, 860)
        self.setStyleSheet(LIGHT_THEME)
        self.editor = None
        self._analysis_context = None
        self._cache_key = self._cache_data = None
        self._current_feature = "pairs"
        self._items_all, self._items_filtered = [], []
        self._page = 0
        self._updating_context = False
        self._loading = False
        self._editor_sizes = [350, 1010]
        self._build_widgets()
        self._context_timer = QtCore.QTimer(self)
        self._context_timer.setInterval(750)
        self._context_timer.timeout.connect(self._check_context)
        self._refresh_objects()

    def _build_widgets(self):
        root = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QHBoxLayout()
        root.addLayout(top)
        self.obj_combo = QtWidgets.QComboBox()
        self.obj_combo.setEditable(True)
        self.obj_combo.setMinimumWidth(220)
        self.obj_combo.currentTextChanged.connect(self._on_object_changed)
        self.state_combo = QtWidgets.QComboBox()
        self.state_combo.currentIndexChanged.connect(self._on_dssr_context_changed)
        self.refresh_obj_btn = DssrUI.button("Refresh objects", self._refresh_objects)
        self.analyze_btn = DssrUI.button(
            "Analyze", lambda: self._load_structure(force=True)
        )
        top.addWidget(QtWidgets.QLabel("Object / selection"))
        top.addWidget(self.obj_combo, 1)
        top.addWidget(QtWidgets.QLabel("State"))
        top.addWidget(self.state_combo)
        top.addWidget(self.refresh_obj_btn)
        top.addWidget(self.analyze_btn)
        self.show_2d_btn = QtWidgets.QPushButton("2D")
        self.show_2d_btn.setCheckable(True)
        self.show_2d_btn.setChecked(True)
        self.show_2d_btn.setToolTip(
            "Show / hide the sequence and 2D panel; keep the edited layout"
        )
        self.show_2d_btn.toggled.connect(self._set_2d_visible)
        top.addWidget(self.show_2d_btn)
        self.dark_btn = QtWidgets.QPushButton("Dark")
        self.dark_btn.setCheckable(True)
        self.dark_btn.setChecked(False)
        self.dark_btn.setToolTip("Toggle PyMOL dark / light mode")
        self.dark_btn.toggled.connect(self._set_dark_mode)
        top.addWidget(self.dark_btn)
        self.settings_btn = QtWidgets.QPushButton("Settings")
        self.settings_btn.setToolTip(
            "Display options, base blocks, and the DSSR executable path"
        )
        self.settings_btn.setCheckable(True)
        top.addWidget(self.settings_btn)

        self.settings_widget = QtWidgets.QWidget()
        settings = QtWidgets.QGridLayout(self.settings_widget)
        settings.setContentsMargins(0, 0, 0, 0)
        settings.setHorizontalSpacing(16)  # Generous separation between columns
        self.exe_edit = QtWidgets.QLineEdit("x3dna-dssr")
        self.exe_edit.textChanged.connect(self._on_dssr_context_changed)
        self.precolor_cb = DssrUI.checkbox("Gray precolor", checked=True)
        self.display_cb = DssrUI.checkbox("Display sticks", checked=False)
        self.zoom_cb = DssrUI.checkbox("Zoom to selection", checked=False)
        self.color_edit = QtWidgets.QLineEdit("auto")
        self.block_file_combo = DssrUI.combo(BLOCK_FEATURES, editable=True)
        self.block_depth_spin = DssrUI.spinbox(0.01, 5.0, 0.5, decimals=True, step=0.05)
        self.make_blocks_btn = DssrUI.button("Make blocks", self._make_blocks_clicked)

        settings.addWidget(QtWidgets.QLabel("DSSR executable"), 0, 0)
        settings.addWidget(self.exe_edit, 0, 1, 1, 5)

        settings.addWidget(self.precolor_cb, 1, 0)
        settings.addWidget(self.display_cb, 1, 1)
        settings.addWidget(self.zoom_cb, 1, 2)

        # Color label + edit paired cleanly
        color_box = QtWidgets.QHBoxLayout()
        color_box.setContentsMargins(0, 0, 0, 0)
        color_box.setSpacing(6)
        color_box.addWidget(QtWidgets.QLabel("Color"))
        color_box.addWidget(self.color_edit, 1)
        settings.addLayout(color_box, 1, 3, 1, 3)

        # Block style: keep the label and dropdown tight together
        block_box = QtWidgets.QHBoxLayout()
        block_box.setContentsMargins(0, 0, 0, 0)
        block_box.setSpacing(6)
        block_box.addWidget(QtWidgets.QLabel("Block style"))
        block_box.addWidget(self.block_file_combo, 1)
        settings.addLayout(block_box, 2, 0, 1, 2)

        # Depth: keep the label and spinbox tight together
        depth_box = QtWidgets.QHBoxLayout()
        depth_box.setContentsMargins(0, 0, 0, 0)
        depth_box.setSpacing(6)
        depth_box.addWidget(QtWidgets.QLabel("Depth"))
        depth_box.addWidget(self.block_depth_spin, 1)
        settings.addLayout(depth_box, 2, 2, 1, 2)

        settings.addWidget(self.make_blocks_btn, 2, 4, 1, 2)

        root.addWidget(self.settings_widget)
        self.settings_widget.hide()
        self.settings_btn.toggled.connect(self.settings_widget.setVisible)

        self.status_label = QtWidgets.QLabel("Load a molecule, then click Analyze.")
        self.status_label.setWordWrap(True)
        root.addWidget(self.status_label)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        root.addWidget(self.splitter, 1)
        sidebar = self.sidebar = QtWidgets.QWidget()
        sidebar.setMinimumWidth(280)
        sidebar.setMaximumWidth(440)
        left = QtWidgets.QVBoxLayout(sidebar)
        left.setContentsMargins(0, 0, 4, 0)
        self.feature_combo = QtWidgets.QComboBox()
        for feature in FEATURE_ORDER:
            self.feature_combo.addItem(FEATURE_LABELS.get(feature, feature), feature)
        self.feature_combo.currentIndexChanged.connect(self._on_feature_changed)
        left.addWidget(self.feature_combo)
        self.filter_edit = QtWidgets.QLineEdit()
        self.filter_edit.setPlaceholderText("Filter (e.g. wc | wobble, -wc -wobble)...")
        self.filter_edit.textChanged.connect(self._on_filter_changed)
        left.addWidget(self.filter_edit)
        self.list_widget = QtWidgets.QListWidget()
        self.list_widget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.list_widget.itemSelectionChanged.connect(
            self._on_selection_changed_preview
        )
        self.list_widget.itemDoubleClicked.connect(self._on_item_double_clicked)
        left.addWidget(self.list_widget, 3)

        paging = QtWidgets.QHBoxLayout()
        self.prev_btn = DssrUI.button("Previous", lambda: self._change_page(-1))
        self.select_all_btn = DssrUI.button(
            "Select all",
            self._select_all_current_feature,
            "Select all items (or all filtered items) in PyMOL",
        )
        self.next_btn = DssrUI.button("Next", lambda: self._change_page(1))
        self.page_label = QtWidgets.QLabel()
        paging.addWidget(self.prev_btn)
        paging.addWidget(self.select_all_btn)
        paging.addWidget(self.page_label, 1)
        paging.addWidget(self.next_btn)
        left.addLayout(paging)

        self.details_box = QtWidgets.QPlainTextEdit()
        self.details_box.setReadOnly(True)
        self.report_box = QtWidgets.QPlainTextEdit()
        self.report_box.setReadOnly(True)
        self.data_tabs = QtWidgets.QTabWidget()
        self.data_tabs.addTab(self.details_box, "Details")
        self.data_tabs.addTab(self.report_box, "Summary")
        self.data_tabs.setCurrentWidget(self.report_box)
        left.addWidget(self.data_tabs, 1)
        self.splitter.addWidget(sidebar)
        self.editor_container = QtWidgets.QWidget()
        self.editor_layout = QtWidgets.QVBoxLayout(self.editor_container)
        self.editor_layout.setContentsMargins(0, 0, 0, 0)
        panel_title = QtWidgets.QHBoxLayout()
        panel_title.addWidget(QtWidgets.QLabel("Sequence · RNA 2D"))
        panel_title.addStretch(1)
        self.minimize_2d_btn = DssrUI.button(
            "−",
            lambda: self.show_2d_btn.setChecked(False),
            "Collapse 2D; keep the layout and undo history",
        )
        self.hide_2d_btn = DssrUI.button(
            "×",
            lambda: self.show_2d_btn.setChecked(False),
            "Hide 2D; reopen with the 2D button above",
        )
        for button in (self.minimize_2d_btn, self.hide_2d_btn):
            button.setFixedWidth(28)
            panel_title.addWidget(button)
        self.editor_layout.addLayout(panel_title)
        self.empty_label = QtWidgets.QLabel(
            "RNA 2D view\n\nLoad a molecule and click Analyze."
        )
        self.empty_label.setAlignment(QtCore.Qt.AlignCenter)
        self.editor_layout.addWidget(self.empty_label)
        self.splitter.addWidget(self.editor_container)
        self.splitter.setStretchFactor(1, 1)
        self.splitter.setSizes([350, 1010])

    def _set_2d_visible(self, visible):
        """Collapse the pane without retiring its analysis or edited coordinates."""
        if not visible:
            self._editor_sizes = self.splitter.sizes()
        self.editor_container.setVisible(visible)
        self.sidebar.setMaximumWidth(440 if visible else 16777215)
        if visible:
            self.splitter.setSizes(self._editor_sizes)
        if self.editor is not None:
            self.editor.set_view_active(visible and not self.isMinimized())
            if visible:
                self.editor.view.setFocus(QtCore.Qt.OtherFocusReason)

    def _set_dark_mode(self, is_dark):
        self.dark_btn.setText("Light" if is_dark else "Dark")
        self.setStyleSheet(DARK_THEME if is_dark else LIGHT_THEME)
        if self.editor is not None:
            self.editor.set_theme(is_dark)

    def changeEvent(self, event):
        super().changeEvent(event)
        if (
            event.type() == QtCore.QEvent.WindowStateChange
            and getattr(self, "editor", None) is not None
        ):
            self.editor.set_view_active(
                self.show_2d_btn.isChecked() and not self.isMinimized()
            )

    def _get_object_text(self):
        return self.obj_combo.currentText().strip() or "all"

    def _get_state_value(self):
        state = self.state_combo.currentData()
        return max(1, int(cmd.get_state())) if state in (None, -1) else int(state)

    def _get_dssr_context(self):
        return (
            self._get_object_text(),
            self.exe_edit.text().strip() or "x3dna-dssr",
            self._get_state_value(),
            int(self.precolor_cb.isChecked()),
        )

    def _context(self):
        selection, exe, state, _precolor = self._get_dssr_context()
        return selection, state, exe

    def _molecule_objects(self):
        """Return a list of valid molecular object names, excluding blocks and selections."""
        objects = cmd.get_object_list()
        valid = []
        for name in objects:
            if (
                name in _DSSR_BLOCK_OBJECTS
                or name in _DSSR_SELECTION_OBJECTS
                or name.startswith("_dssr_2d_")
            ):
                continue
            try:
                # Ensure the entry is a real molecule and contains atoms
                obj_type = cmd.get_type(name)
                if obj_type == "object:molecule" and cmd.count_atoms(name) > 0:
                    valid.append(name)
            except Exception:
                # Fallback check if get_type is unavailable
                try:
                    if cmd.count_atoms(name) > 0:
                        valid.append(name)
                except Exception:
                    pass
        return valid

    def _update_state_combo(self, wanted=None):
        selected = self.state_combo.currentData() if wanted is None else wanted
        try:
            count = max(1, int(cmd.count_states(self._get_object_text())))
        except Exception:
            count = 1
        self.state_combo.blockSignals(True)
        self.state_combo.clear()
        self.state_combo.addItem("Current", -1)
        for state in range(1, count + 1):
            self.state_combo.addItem(str(state), state)
        index = self.state_combo.findData(selected)
        self.state_combo.setCurrentIndex(max(0, index))
        self.state_combo.blockSignals(False)

    def _refresh_objects(self):
        objects = self._molecule_objects()
        previous = self.obj_combo.currentText().strip()
        try:
            keep_previous = bool(previous and cmd.count_atoms(previous) > 0)
        except Exception:
            keep_previous = False
        self._updating_context = True
        try:
            self.obj_combo.clear()
            self.obj_combo.addItems(objects)
            if keep_previous and objects:
                self.obj_combo.setEditText(previous)
            elif objects:
                self.obj_combo.setCurrentIndex(0)
            else:
                self.obj_combo.setEditText("")
            self._update_state_combo()
        finally:
            self._updating_context = False
        if not objects:
            self._clear_analysis(
                "No molecule loaded. Load a PDB/CIF file, then click Analyze."
            )
        elif self._analysis_context != self._context():
            self._clear_analysis("Choose an object / selection, then click Analyze.")

    def _on_object_changed(self, *_args):
        if not self._updating_context:
            self._update_state_combo()
            self._on_dssr_context_changed()

    def _on_dssr_context_changed(self, *_args):
        if not self._updating_context:
            self._clear_analysis(
                "Analysis context changed. Click Analyze to update the structure."
            )

    def _invalidate_dssr_cache(self):
        DssrUtils._invalidate_cache()

    def _get_dssr_data(self, selection, state, exe, precolor_on):
        if not self._molecule_objects() or cmd.count_atoms(selection, state=state) <= 0:
            raise CmdException("No atoms in the requested object / selection.")
        return DssrUtils._cached_selection_json(
            selection, state, exe, precolor=bool(precolor_on)
        )

    def _dispose_editor(self):
        if self.editor is not None:
            editor, self.editor = self.editor, None
            editor.shutdown()
            self.editor_layout.removeWidget(editor)
            editor.hide()
            editor.deleteLater()

    def _clear_analysis(self, message):
        self._dispose_editor()
        self._analysis_context = None
        self._invalidate_dssr_cache()
        self._items_all, self._items_filtered = [], []
        self._page = 0
        self.details_box.clear()
        self.report_box.clear()
        self._update_feature_counts(None)
        self._render_list()
        self.empty_label.setText(message)
        self.empty_label.show()
        self.status_label.setText(message)

        # Clear temporary selection indicators in PyMOL
        try:
            cmd.delete("sele")
            cmd.delete("indicate")
            cmd.refresh()
        except Exception:
            pass

    def _check_context(self):
        if self._loading or self._analysis_context is None or not self.isVisible():
            return
        try:
            selection, state, _exe = self._context()
            if (
                not self._molecule_objects()
                or cmd.count_atoms(selection, state=state) <= 0
            ):
                self._clear_analysis(
                    "The analyzed structure is no longer loaded. Click Analyze after loading it."
                )
            elif self._analysis_context != self._context():
                self._on_dssr_context_changed()
        except Exception as error:
            self._clear_analysis("Structure context is unavailable: %s" % error)

    def _require_analysis(self):
        self._check_context()
        data = _DSSR_DATA_CACHE.get("data")
        if self._analysis_context != self._context() or data is None:
            raise CmdException("Click Analyze for the current object and state first.")
        return data

    def _big_object_warning(self, sel):
        """Return a warning string if the selection exceeds 20,000 atoms."""
        thresh = 20000
        try:
            n_atoms = int(cmd.count_atoms(sel))
        except Exception:
            n_atoms = 0
        if n_atoms >= thresh:
            return (
                "Warning: Large selection (%d atoms). DSSR analysis may be slow. "
                "Consider selecting a specific chain.\n" % n_atoms
            )
        return ""

    def _load_structure(self, force=True):
        if self._loading:
            return self.editor
        selection, exe, state, precolor = self._get_dssr_context()
        if (
            not force
            and self.editor is not None
            and self._analysis_context == (selection, state, exe)
        ):
            return self.editor

        self._loading = True
        self.analyze_btn.setEnabled(False)
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        warn = self._big_object_warning(selection)
        self.status_label.setText(
            warn + ("Analyzing %s, state %d..." % (selection, state))
        )
        QtWidgets.QApplication.processEvents()

        try:
            if force:
                self._dispose_editor()
                self._invalidate_dssr_cache()
            data = self._get_dssr_data(selection, state, exe, precolor)
            return self.show_analysis(data, selection, state, exe, force=force)
        except Exception as error:
            self._clear_analysis("Analysis error: %s" % error)
            return None
        finally:
            self._loading = False
            self.analyze_btn.setEnabled(True)
            QtWidgets.QApplication.restoreOverrideCursor()

    def show_analysis(
        self,
        data,
        selection,
        state,
        exe,
        algorithm="standard",
        number_every=10,
        show_tertiary=0,
        title="",
        force=False,
    ):
        context = (str(selection), int(state), str(exe))
        if not force and self.editor is not None and context == self._analysis_context:
            if not self._loading:
                self._update_state_combo(wanted=int(state))
            requested = str(algorithm or "standard").strip().lower()
            requested = requested if requested in LAYOUT_CHOICES else "standard"
            if self.editor.layout_combo.currentText() != requested:
                self.editor.layout_combo.setCurrentText(requested)
            self.editor.number_spin.setValue(max(0, int(number_every)))
            self.editor.tertiary_cb.setChecked(bool(int(show_tertiary)))
            if title:
                self.editor.model.title = str(title)
            return self.editor
        model = Dssr2DModel.from_dssr(
            data, title=title or "%s — state %d" % (selection, state)
        )
        self._dispose_editor()
        keep_current = self._loading and self.state_combo.currentData() == -1
        self._updating_context = True
        try:
            self.obj_combo.setEditText(str(selection))
            self.exe_edit.setText(str(exe))
            self._update_state_combo(wanted=-1 if keep_current else int(state))
        finally:
            self._updating_context = False
        self._analysis_context = self._cache_key = context
        self._cache_data = data
        self.editor = Dssr2DEditor(
            model, selection, algorithm, number_every, show_tertiary, parent=self
        )
        self.editor.pymol_state = int(state)
        self.editor.set_theme(self.dark_btn.isChecked())
        self.editor_layout.addWidget(self.editor)
        self.empty_label.hide()
        self.editor.show()
        self.editor.set_view_active(
            self.show_2d_btn.isChecked() and not self.isMinimized()
        )
        if self.editor.isVisible():
            self.editor.view.setFocus(QtCore.Qt.OtherFocusReason)
        self._update_feature_counts(data)
        self._refresh_list()
        self.report_box.setPlainText(DssrParser._format_rna_summary_text(data))
        self.data_tabs.setCurrentWidget(self.report_box)
        self.status_label.setText(
            "%s | state %d | %s" % (selection, state, model.summary())
        )
        QtCore.QTimer.singleShot(0, self.editor.fit_scene)
        return self.editor

    def _update_feature_counts(self, data):
        current_feat = self.feature_combo.currentData()
        self.feature_combo.blockSignals(True)
        self.feature_combo.clear()

        if data is not None:
            for feature in FEATURE_ORDER:
                count = (
                    DssrParser._count_pseudoknot_layers(data)
                    if feature == "pseudoknot"
                    else len(DssrParser.feature_entries(data, feature))
                )
                if count > 0:
                    label = "%s (%d)" % (FEATURE_LABELS.get(feature, feature), count)
                    self.feature_combo.addItem(label, feature)

        idx = self.feature_combo.findData(current_feat)
        if idx >= 0:
            self.feature_combo.setCurrentIndex(idx)
        elif self.feature_combo.count() > 0:
            self.feature_combo.setCurrentIndex(0)

        has_items = self.feature_combo.count() > 0
        self.feature_combo.setEnabled(has_items)
        self.make_blocks_btn.setEnabled(has_items)
        self.feature_combo.blockSignals(False)

        self._current_feature = str(self.feature_combo.currentData() or "pairs")

    def _on_feature_changed(self, *_args):
        self._current_feature = str(self.feature_combo.currentData() or "pairs")
        self.details_box.clear()
        self._refresh_list()

    def _refresh_list(self):
        self._items_all = []
        data = _DSSR_DATA_CACHE.get("data")
        if data is not None:
            if self._current_feature == "pseudoknot":
                layers = DssrParser.parse_dotbracket_pseudoknots(
                    DssrParser._extract_dotbracket(data)
                )
                self._items_all = [
                    (i, "%d: layer %s, %d pairs" % (i, key, len(layers[key])))
                    for i, key in enumerate(sorted(layers), 1)
                ]
            else:
                self._items_all = [
                    (i, self._entry_label(self._current_feature, entry, i))
                    for i, entry in enumerate(
                        DssrParser.feature_entries(data, self._current_feature), 1
                    )
                ]
        self._page = 0
        self._render_list()

    @staticmethod
    def _entry_label(feature, entry, index):
        """Compact display labels only; DSSR identifiers and item indices stay intact."""

        def short_id(match):
            identifier = match.group(0)
            try:
                chain, resi = DssrParser.parse_nt_id(identifier)
                parts = identifier.split("|")
                atom = "/" + parts[5] if len(parts) > 5 and parts[5] else ""
                return "%s:%s%s%s" % (chain or "–", parts[3], resi, atom)
            except (CmdException, IndexError):
                return identifier

        text = DssrParser._preview_entry(feature, entry, index)
        return re.sub(r"[^\s,;()]+\|[^\s,;()]*", short_id, text)

    def _on_filter_changed(self, *_args):
        self._page = 0
        self._render_list()

    def _change_page(self, delta):
        self._page = max(0, self._page + delta)
        self._render_list()

    def _render_list(self):
        query = self.filter_edit.text().strip()
        self._items_filtered = [
            item
            for item in self._items_all
            if DssrUtils.matches_boolean_query(item[1], query)
        ]
        count = len(self._items_filtered)
        pages = max(1, (count + self.PAGE_SIZE - 1) // self.PAGE_SIZE)
        self._page = min(self._page, pages - 1)
        self.list_widget.clear()
        start = self._page * self.PAGE_SIZE
        for index, text in self._items_filtered[start : start + self.PAGE_SIZE]:
            item = QtWidgets.QListWidgetItem(text)
            item.setData(QtCore.Qt.UserRole, index)
            self.list_widget.addItem(item)
        self.page_label.setText("%d items · %d/%d" % (count, self._page + 1, pages))
        self.prev_btn.setEnabled(self._page > 0)
        self.next_btn.setEnabled(self._page + 1 < pages)

    def _entry_for_feature_index(self, data, feature, index):
        if feature == "pseudoknot":
            layers = DssrParser.parse_dotbracket_pseudoknots(
                DssrParser._extract_dotbracket(data)
            )
            key = sorted(layers)[index - 1]
            return {"layer": key, "pair_count": len(layers[key]), "pairs": layers[key]}
        return DssrParser.feature_entries(data, feature)[index - 1]

    def _show_item_details(self, data, feature, index):
        entry = self._entry_for_feature_index(data, feature, index)
        lines = ["%s #%d" % (FEATURE_LABELS.get(feature, feature), index)]
        for key, value in entry.items():
            if isinstance(value, (list, dict)):
                value = "%d entries" % len(value)
            lines.append("%s: %s" % (key, str(value)[:200]))
        self.details_box.setPlainText("\n".join(lines))
        self.data_tabs.setCurrentWidget(self.details_box)

    def _on_selection_changed_preview(self):
        items = self.list_widget.selectedItems()

        # Clear the temporary PyMOL selection and 2D highlights if nothing is selected
        if not items:
            try:
                cmd.delete("sele")
                cmd.refresh()
            except Exception:
                pass
            self.details_box.clear()
            if self.editor is not None:
                self.editor.clear_base_selection()
            return

        try:
            data = self._require_analysis()
            selection = self._analysis_context[0]
            cores = []

            for item in items:
                index = item.data(QtCore.Qt.UserRole)
                if index is not None:
                    core = DssrParser._build_residue_sel_from_dssr(
                        data, self._current_feature, int(index)
                    )
                    if core:
                        cores.append(core)

            if cores:
                sel_str = " or ".join("(%s)" % c for c in cores)

                cmd.select("sele", "byres ((%s) and (%s))" % (selection, sel_str))
                cmd.enable("sele")
                cmd.refresh()

                if len(items) == 1:
                    index = items[0].data(QtCore.Qt.UserRole)
                    self._show_item_details(data, self._current_feature, int(index))
                else:
                    self.details_box.setPlainText(
                        "%d items selected for preview." % len(items)
                    )
                    self.data_tabs.setCurrentWidget(self.details_box)

                if self.editor is not None:
                    all_residues = set()
                    cmd.iterate(
                        "sele",
                        "_dssr_res.add((chain, resi))",
                        space={"_dssr_res": all_residues},
                    )
                    matching = {
                        i
                        for i, nt in enumerate(self.editor.model.nts)
                        if (
                            str(nt.get("chain", "")).strip(),
                            str(nt.get("resi", "")).strip(),
                        )
                        in all_residues
                    }

                    blocked = self.editor.scene.blockSignals(True)
                    was_rebuilding = self.editor._rebuilding
                    self.editor._rebuilding = True
                    try:
                        for node in self.editor.nodes:
                            node.setSelected(node.nt_index in matching)
                    finally:
                        self.editor._rebuilding = was_rebuilding
                        self.editor.scene.blockSignals(blocked)

                    self.editor._sync_sequence_selection()
                    self.editor._last_pymol_signature = tuple()
                    self.editor._sync_timer.stop()
                    self.editor._sync_pending = False

        except Exception as error:
            self.status_label.setText("Selection error: %s" % error)

    def _on_item_double_clicked(self, item):
        data = item.data(QtCore.Qt.UserRole)
        if data is None:
            return
        try:
            idx = int(data)
        except Exception:
            return

        sel = self._get_object_text()
        feat = self._current_feature
        exe = self.exe_edit.text().strip() or "x3dna-dssr"
        st = self._get_state_value()

        # Build standard selection name (e.g., junctions1, stems2, uturns1)
        nm = "%s%d" % (feat.lower(), idx)

        col = self.color_edit.text().strip() or "auto"
        precolor_on = 1 if self.precolor_cb.isChecked() else 0
        display_on = 1 if self.display_cb.isChecked() else 0
        zoom_on = 1 if self.zoom_cb.isChecked() else 0
        showinfo_on = 0
        radius = 0.25

        try:
            DssrCmd._dssr(
                sel=sel,
                f=feat,
                i=idx,
                n=nm,
                q=0,
                si=showinfo_on,
                st=st,
                exe=exe,
                color=col,
                display=display_on,
                stick_radius=radius,
                do_zoom=zoom_on,
                pc=precolor_on,
            )

            # Drop temporary (sele) without deselecting the named object
            try:
                cmd.delete("sele")
                cmd.delete("indicate")
            except Exception:
                pass

            # Force state transition so PyMOL paints selection markers
            cmd.disable(nm)
            cmd.enable(nm)

            # Update the 2D layout canvas nodes silently
            if self.editor is not None:
                all_residues = set()
                cmd.iterate(
                    nm,
                    "_dssr_res.add((chain, resi))",
                    space={"_dssr_res": all_residues},
                )
                matching = {
                    i
                    for i, nt in enumerate(self.editor.model.nts)
                    if (
                        str(nt.get("chain", "")).strip(),
                        str(nt.get("resi", "")).strip(),
                    )
                    in all_residues
                }

                blocked = self.editor.scene.blockSignals(True)
                was_rebuilding = self.editor._rebuilding
                self.editor._rebuilding = True
                try:
                    for node in self.editor.nodes:
                        node.setSelected(node.nt_index in matching)
                finally:
                    self.editor._rebuilding = was_rebuilding
                    self.editor.scene.blockSignals(blocked)

                self.editor._sync_sequence_selection()
                self.editor._last_pymol_signature = tuple()

        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, "DSSR GUI error", str(e))
            except Exception:
                pass
            try:
                print("dssr_gui select error: %s" % str(e))
            except Exception:
                pass

    def _select_all_current_feature(self):
        try:
            data = self._require_analysis()
            selection, state, _exe = self._analysis_context
            feature = self._current_feature

            has_filter = bool(self.filter_edit.text().strip())

            if has_filter and self._items_filtered:
                indices = [item[0] for item in self._items_filtered]
            elif feature == "pseudoknot":
                layers = DssrParser.parse_dotbracket_pseudoknots(
                    DssrParser._extract_dotbracket(data)
                )
                indices = list(range(1, len(layers) + 1))
            else:
                entries = DssrParser.feature_entries(data, feature)
                indices = list(range(1, len(entries) + 1))

            if not indices:
                self.status_label.setText("No %s items found to select." % feature)
                return

            self.list_widget.blockSignals(True)
            self.list_widget.selectAll()
            self.list_widget.blockSignals(False)

            parts = [
                DssrParser._build_residue_sel_from_dssr(data, feature, idx)
                for idx in indices
            ]
            sel_str = " or ".join("(%s)" % part for part in parts if part)
            if not sel_str:
                self.status_label.setText("Could not build selection for %s." % feature)
                return

            name = "%s_%s" % (feature.lower(), "filtered" if has_filter else "all")
            DssrCmd._create_feature_selection(name, selection, sel_str, quiet=0)

            col = self.color_edit.text().strip() or "auto"
            user_color = DssrUtils._resolve_color_spec(col)
            cmd.color(user_color if user_color else "pink", name)

            if self.display_cb.isChecked():
                DssrCmd._display_feature_selection(
                    name, display=1, stick_radius=0.25, do_zoom=0
                )

            cmd.select("sele", name)
            cmd.enable("sele")
            cmd.refresh()

            if self.zoom_cb.isChecked():
                cmd.zoom(name)

            self.status_label.setText(
                "Created selection '%s' with %d %s items."
                % (name, len(indices), feature)
            )

            if self.editor is not None:
                all_residues = set()
                cmd.iterate(
                    name,
                    "_dssr_res.add((chain, resi))",
                    space={"_dssr_res": all_residues},
                )
                matching = {
                    idx
                    for idx, nt in enumerate(self.editor.model.nts)
                    if (
                        str(nt.get("chain", "")).strip(),
                        str(nt.get("resi", "")).strip(),
                    )
                    in all_residues
                }

                blocked = self.editor.scene.blockSignals(True)
                was_rebuilding = self.editor._rebuilding
                self.editor._rebuilding = True
                try:
                    for node in self.editor.nodes:
                        node.setSelected(node.nt_index in matching)
                finally:
                    self.editor._rebuilding = was_rebuilding
                    self.editor.scene.blockSignals(blocked)

                self.editor._sync_sequence_selection()
                self.editor._last_pymol_signature = tuple(sorted(all_residues))
                self.editor._update_editor_status(
                    "%d %s selected" % (len(indices), feature)
                )

        except Exception as error:
            self.status_label.setText("Select all error: %s" % error)

    def _make_blocks_clicked(self):
        try:
            data = self._require_analysis()
            selection, state, exe = self._analysis_context
            indices = [
                item.data(QtCore.Qt.UserRole)
                for item in self.list_widget.selectedItems()
            ]
            parts = [
                DssrParser._build_residue_sel_from_dssr(
                    data, self._current_feature, int(index)
                )
                for index in indices
                if index is not None
            ]

            if parts:
                core = " or ".join("(%s)" % part for part in parts)
            else:
                feature_all_name = "%s_all" % self._current_feature.lower()
                if (
                    feature_all_name in cmd.get_names("selections")
                    and cmd.count_atoms(feature_all_name) > 0
                ):
                    core = feature_all_name
                elif (
                    "sele" in cmd.get_names("selections")
                    and cmd.count_atoms("sele") > 0
                ):
                    core = "sele"
                elif self.editor is not None and any(
                    n.isSelected() for n in self.editor.nodes
                ):
                    sig = self.editor._node_residue_signature()
                    res_parts = DssrParser._compact_sel_from_residues(
                        {(c, r) for c, r in sig if c}
                    )
                    core = res_parts if res_parts else "all"
                else:
                    raise CmdException("Select a feature or some bases first.")

            scope = "byres ((%s) and (%s))" % (selection, core)
            if cmd.count_atoms(scope, state=state) <= 0:
                raise CmdException("Select a feature or some bases first.")

            name = DssrCmd._unused_name("dssr_blocks")
            DssrCmd.dssr_block(
                selection=scope,
                state=state,
                block_file=self.block_file_combo.currentText().strip() or "face",
                block_depth=float(self.block_depth_spin.value()),
                name=name,
                exe=exe,
                quiet=1,
            )
            if self.zoom_cb.isChecked():
                cmd.zoom(name)
            self.status_label.setText("Created %s for the selected bases." % name)
        except Exception as error:
            self.status_label.setText("Blocks error: %s" % error)

    def showEvent(self, event):
        super().showEvent(event)
        self._context_timer.start()

    def closeEvent(self, event):
        self._context_timer.stop()
        self._clear_analysis("Click Analyze to load the structure again.")
        super().closeEvent(event)

    def reject(self):
        self._context_timer.stop()
        self._clear_analysis("Click Analyze to load the structure again.")
        super().reject()

    @staticmethod
    def dssr_gui():
        global _DSSR_GUI_DIALOG

        has_structure = bool(
            [
                name
                for name in cmd.get_object_list()
                if name not in _DSSR_BLOCK_OBJECTS
                and not name.startswith("_dssr_2d_")
                and cmd.count_atoms(name) > 0
            ]
        )

        if _DSSR_GUI_DIALOG is None:
            _DSSR_GUI_DIALOG = DssrGuiDialog()
        host = _DSSR_GUI_DIALOG
        host.show()
        host.raise_()
        host.activateWindow()
        host._check_context()

        if host.editor is None:
            host._refresh_objects()
            if has_structure:
                host._load_structure(force=False)

        if not has_structure:
            msg = "No structure loaded. Please load a PDB/CIF file before running DSSR-PyMOL"
            print(msg)
            host.status_label.setText(msg)
            QtWidgets.QMessageBox.warning(host, "No Structure Loaded", msg)

        if host.editor is not None and host.editor.isVisible():
            host.editor.view.setFocus(QtCore.Qt.OtherFocusReason)
        return host


# RNA data model
class Dssr2DModel:
    """Normalized sequence, dot-bracket, residue, and pairing information."""

    _OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
    _CLOSE_TO_OPEN = {v: k for k, v in _OPEN_TO_CLOSE.items()}
    _LAYER_BY_OPEN = {"(": 0, "[": 1, "{": 2, "<": 3}

    def __init__(self):
        self.sequence = ""
        self.structure = ""
        self.raw_sequence = ""
        self.raw_structure = ""
        self.nts = []
        self.chain_breaks = set()
        self.secondary_pairs = []
        self.tertiary_pairs = []
        self.warnings = []
        self.title = "RNA secondary structure"

    @staticmethod
    def _clean_dbn_text(value, keep_ampersand=True):
        if value is None:
            return ""
        text = str(value)
        out = []
        for ch in text:
            if ch.isspace():
                continue
            if ch == "&" and keep_ampersand:
                out.append(ch)
                continue
            out.append(ch)
        return "".join(out)

    @staticmethod
    def _find_dbn_record(dssr_data):
        dbn = dssr_data.get("dbn", None)
        if isinstance(dbn, dict):
            all_chains = dbn.get("all_chains")
            if isinstance(all_chains, dict):
                return all_chains
            if any(k in dbn for k in ("sstr", "bseq", "seq", "sequence")):
                return dbn
            for value in dbn.values():
                if isinstance(value, dict) and any(
                    k in value for k in ("sstr", "bseq", "seq", "sequence")
                ):
                    return value
        if isinstance(dbn, str):
            return {"sstr": dbn}
        return {}

    @staticmethod
    def _base_from_nt_entry(nt):
        if not isinstance(nt, dict):
            return "N"
        for key in (
            "nt_code",
            "base",
            "bseq",
            "one_letter",
            "nt_name",
            "resname",
        ):
            value = nt.get(key)
            if value is None:
                continue
            text = str(value).strip()
            if not text:
                continue
            if len(text) == 1:
                return text
            upper = text.upper()
            for base in ("A", "C", "G", "U", "T", "I"):
                if upper == base or upper.endswith(base):
                    return base
            return text[0]
        return "N"

    @staticmethod
    def _safe_parse_nt_id(nt_id):
        try:
            chain, resi = DssrParser.parse_nt_id(nt_id)
            return str(chain), str(resi)
        except Exception:
            return "", ""

    @classmethod
    def from_dssr(cls, dssr_data, title="RNA secondary structure"):
        if not isinstance(dssr_data, dict):
            raise CmdException("DSSR JSON data must be a dictionary")

        model = cls()
        model.title = str(title or "RNA secondary structure")

        record = cls._find_dbn_record(dssr_data)
        seq_raw = ""
        for key in ("bseq", "seq", "sequence"):
            if record.get(key) is not None:
                seq_raw = cls._clean_dbn_text(record.get(key), keep_ampersand=True)
                if seq_raw:
                    break

        sstr_raw = ""
        for key in ("sstr", "structure", "dbn"):
            if record.get(key) is not None:
                sstr_raw = cls._clean_dbn_text(record.get(key), keep_ampersand=True)
                if sstr_raw:
                    break

        if not sstr_raw:
            try:
                sstr_raw = cls._clean_dbn_text(
                    DssrParser._extract_dotbracket(dssr_data),
                    keep_ampersand=True,
                )
            except Exception:
                sstr_raw = ""

        nts_json = dssr_data.get("nts", [])
        if not isinstance(nts_json, list):
            nts_json = []

        seq_no_breaks = seq_raw.replace("&", "")
        sstr_no_breaks = sstr_raw.replace("&", "")
        nts_count = len(nts_json)

        if not seq_no_breaks and nts_count:
            seq_no_breaks = "".join(cls._base_from_nt_entry(nt) for nt in nts_json)

        target_len = 0
        if sstr_no_breaks:
            target_len = len(sstr_no_breaks)
        elif seq_no_breaks:
            target_len = len(seq_no_breaks)
        elif nts_count:
            target_len = nts_count

        if target_len <= 0:
            raise CmdException(
                "DSSR output contains no usable sequence or dot-bracket structure"
            )

        if nts_count and target_len != nts_count:
            # DSSR's nts array is the most useful source for PyMOL residue mapping.
            # Prefer it when the DBN record omitted separators or modified bases.
            if not sstr_no_breaks or abs(target_len - nts_count) <= 2:
                target_len = nts_count

        if len(seq_no_breaks) < target_len:
            derived = "".join(cls._base_from_nt_entry(nt) for nt in nts_json)
            if len(derived) >= target_len:
                seq_no_breaks = derived[:target_len]
            else:
                seq_no_breaks = (seq_no_breaks + derived + ("N" * target_len))[
                    :target_len
                ]
        elif len(seq_no_breaks) > target_len:
            seq_no_breaks = seq_no_breaks[:target_len]

        if len(sstr_no_breaks) < target_len:
            sstr_no_breaks = sstr_no_breaks + ("." * (target_len - len(sstr_no_breaks)))
        elif len(sstr_no_breaks) > target_len:
            sstr_no_breaks = sstr_no_breaks[:target_len]

        model.raw_sequence = seq_raw
        model.raw_structure = sstr_raw
        model.sequence = seq_no_breaks
        model.structure = sstr_no_breaks

        # Chain breaks explicitly encoded in sequence or structure.
        for raw in (seq_raw, sstr_raw):
            nt_pos = 0
            for ch in raw:
                if ch == "&":
                    if nt_pos > 0:
                        model.chain_breaks.add(nt_pos - 1)
                else:
                    nt_pos += 1

        # Normalize nucleotide metadata to exactly target_len entries.
        previous_chain = None
        for i in range(target_len):
            src = (
                nts_json[i]
                if i < len(nts_json) and isinstance(nts_json[i], dict)
                else {}
            )
            nt_id = str(src.get("nt_id", ""))
            chain, resi = cls._safe_parse_nt_id(nt_id)
            base = (
                model.sequence[i]
                if i < len(model.sequence)
                else cls._base_from_nt_entry(src)
            )
            nt = {
                "index": i,
                "number": i + 1,
                "base": str(base),
                "nt_id": nt_id,
                "chain": chain,
                "resi": resi,
                "name": str(src.get("nt_name", src.get("resname", ""))),
                "source": src,
            }
            model.nts.append(nt)
            if previous_chain is not None and chain and chain != previous_chain:
                model.chain_breaks.add(i - 1)
            if chain:
                previous_chain = chain

        model.secondary_pairs = model._parse_dotbracket_pairs(model.structure)
        model._merge_dssr_pairs(dssr_data)
        return model

    def _parse_dotbracket_pairs(self, structure):
        stacks = {}
        pairs = []
        letter_layers = {}
        next_layer = 4

        for idx, char in enumerate(str(structure)):
            if char in self._OPEN_TO_CLOSE:
                stacks.setdefault(char, []).append(idx)
                continue

            if char in self._CLOSE_TO_OPEN:
                opener = self._CLOSE_TO_OPEN[char]
                stack = stacks.setdefault(opener, [])
                if not stack:
                    self.warnings.append(
                        "Unmatched closing bracket %s at nucleotide %d"
                        % (char, idx + 1)
                    )
                    continue
                i = stack.pop()
                layer = self._LAYER_BY_OPEN.get(opener, 0)
                pairs.append(
                    {
                        "i": i,
                        "j": idx,
                        "layer": layer,
                        "symbol": opener + char,
                        "lw": "",
                        "kind": "secondary",
                    }
                )
                continue

            # Extended Vienna notation: A..a, B..b, ... represent additional
            # pseudoknot layers.
            if "A" <= char <= "Z":
                stacks.setdefault(char, []).append(idx)
                if char not in letter_layers:
                    letter_layers[char] = next_layer
                    next_layer += 1
                continue

            if "a" <= char <= "z":
                opener = char.upper()
                stack = stacks.setdefault(opener, [])
                if not stack:
                    self.warnings.append(
                        "Unmatched pseudoknot symbol %s at nucleotide %d"
                        % (char, idx + 1)
                    )
                    continue
                i = stack.pop()
                layer = letter_layers.setdefault(opener, next_layer)
                if layer == next_layer:
                    next_layer += 1
                pairs.append(
                    {
                        "i": i,
                        "j": idx,
                        "layer": layer,
                        "symbol": opener + char,
                        "lw": "",
                        "kind": "secondary",
                    }
                )

        for opener, stack in stacks.items():
            for idx in stack:
                self.warnings.append(
                    "Unmatched opening symbol %s at nucleotide %d" % (opener, idx + 1)
                )

        pairs.sort(key=lambda p: (p["i"], p["j"]))
        return pairs

    def _merge_dssr_pairs(self, dssr_data):
        by_exact_id = {}
        by_chain_resi = {}
        for nt in self.nts:
            if nt["nt_id"]:
                by_exact_id[nt["nt_id"]] = nt["index"]
            if nt["chain"] or nt["resi"]:
                by_chain_resi[(nt["chain"], nt["resi"])] = nt["index"]

        secondary_by_key = {}
        for pair in self.secondary_pairs:
            key = tuple(sorted((int(pair["i"]), int(pair["j"]))))
            secondary_by_key[key] = pair

        tertiary_by_key = {}
        entries = dssr_data.get("pairs", [])
        if not isinstance(entries, list):
            entries = []

        for entry in entries:
            if not isinstance(entry, dict):
                continue
            nt1 = str(entry.get("nt1", ""))
            nt2 = str(entry.get("nt2", ""))
            i = by_exact_id.get(nt1)
            j = by_exact_id.get(nt2)
            if i is None and nt1:
                i = by_chain_resi.get(self._safe_parse_nt_id(nt1))
            if j is None and nt2:
                j = by_chain_resi.get(self._safe_parse_nt_id(nt2))
            if i is None or j is None or i == j:
                continue

            key = tuple(sorted((int(i), int(j))))
            lw = str(entry.get("LW", entry.get("bp", "")))
            if key in secondary_by_key:
                secondary_by_key[key]["lw"] = lw
                secondary_by_key[key]["dssr"] = entry
                continue

            if key not in tertiary_by_key:
                tertiary_by_key[key] = {
                    "i": key[0],
                    "j": key[1],
                    "layer": -1,
                    "symbol": "",
                    "lw": lw,
                    "kind": "tertiary",
                    "dssr": entry,
                }

        self.tertiary_pairs = sorted(
            tertiary_by_key.values(), key=lambda p: (p["i"], p["j"])
        )

    def planar_secondary_pairs(self):
        """Return a deterministic non-crossing scaffold for layout.

        Standard parenthesis pairs are preferred, then additional dot-bracket
        layers are admitted only when they neither reuse an endpoint nor cross
        a pair already in the scaffold. All excluded pairs are still rendered
        as pseudoknot/auxiliary edges.
        """
        candidates = sorted(
            self.secondary_pairs,
            key=lambda p: (
                0 if int(p.get("layer", 0)) == 0 else 1,
                int(p.get("layer", 0)),
                -(int(p.get("j", 0)) - int(p.get("i", 0))),
                int(p.get("i", 0)),
            ),
        )
        scaffold = []
        used = set()
        for pair in candidates:
            i = int(pair.get("i", -1))
            j = int(pair.get("j", -1))
            if i > j:
                i, j = j, i
            if i < 0 or j < 0 or i == j or i in used or j in used:
                continue
            crossing = False
            for accepted in scaffold:
                a = int(accepted["i"])
                b = int(accepted["j"])
                if a > b:
                    a, b = b, a
                if (a < i < b < j) or (i < a < j < b):
                    crossing = True
                    break
            if crossing:
                continue
            scaffold.append(pair)
            used.add(i)
            used.add(j)
        scaffold.sort(key=lambda p: (int(p["i"]), int(p["j"])))
        return scaffold

    def chain_count(self):
        chains = []
        seen = set()
        for nt in self.nts:
            chain = nt.get("chain", "")
            if chain and chain not in seen:
                seen.add(chain)
                chains.append(chain)
        if chains:
            return len(chains)
        return len(self.chain_breaks) + 1

    def summary(self):
        return (
            "%d nt | %d chain(s) | %d secondary pair(s) | %d additional DSSR pair(s)"
            % (
                len(self.nts),
                self.chain_count(),
                len(self.secondary_pairs),
                len(self.tertiary_pairs),
            )
        )


# ---------------------------------------------------------------------------
# Bundled NAView geometry engine (1-based internal indexing, planar scaffold)
# ---------------------------------------------------------------------------
# Python adaptation of ViennaRNA/fornac src/naview/naview.js
# Fornac authors: Peter Kerpedjiev, Stefan Hammer, and Ronny Lorenz
# Licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
# ---------------------------------------------------------------------------


class DssrNaviewRegion:
    __slots__ = ("start1", "end1", "start2", "end2")

    def __init__(self):
        self.start1 = 0
        self.end1 = 0
        self.start2 = 0
        self.end2 = 0


class DssrNaviewBase:
    __slots__ = ("mate", "x", "y", "extracted", "region")

    def __init__(self):
        self.mate = 0
        self.x = 9999.0
        self.y = 9999.0
        self.extracted = False
        self.region = None


class DssrNaviewConnection:
    __slots__ = (
        "loop",
        "region",
        "start",
        "end",
        "xrad",
        "yrad",
        "angle",
        "extruded",
    )

    def __init__(self):
        self.loop = None
        self.region = None
        self.start = 0
        self.end = 0
        self.xrad = 0.0
        self.yrad = 0.0
        self.angle = 0.0
        self.extruded = False


class DssrNaviewLoop:
    __slots__ = (
        "connections",
        "depth",
        "mark",
        "x",
        "y",
        "radius",
    )

    def __init__(self):
        self.connections = []
        self.depth = 0
        self.mark = False
        self.x = 0.0
        self.y = 0.0
        self.radius = 0.0

    @property
    def nconnection(self):
        return len(self.connections)


class DssrNaview:
    """Dependency-free Python NAView coordinates for a planar RNA scaffold."""

    ANUM = 9999.0
    MAXITER = 500
    HELIX_FACTOR = 0.6
    BACKBONE_DISTANCE = 27.0

    def __init__(self):
        self.nbase = 0
        self.nregion = 0
        self.loop_count = 0
        self.root = None
        self.bases = []
        self.regions = []
        self.loops = []
        self.lencut = 0.8
        self.RADIUS_REDUCTION_FACTOR = 1.4
        self.angleinc = 0.0
        self._h = 0.0

    def coordinates(self, pair_table):
        if not pair_table or int(pair_table[0]) <= 0:
            return []

        self.nbase = int(pair_table[0])
        if len(pair_table) != self.nbase + 1:
            raise ValueError(
                "NAView pair table length %d does not match n=%d"
                % (len(pair_table), self.nbase)
            )

        self.bases = [DssrNaviewBase() for _ in range(self.nbase + 1)]
        self.regions = [DssrNaviewRegion() for _ in range(self.nbase + 1)]
        self._read_in_bases(pair_table)
        self._find_regions()

        self.loops = [DssrNaviewLoop() for _ in range(self.nbase + 1)]
        self.loop_count = 0
        self._construct_loop(0)
        self._find_central_loop()
        if self.root is None:
            raise RuntimeError("NAView could not identify a central loop")
        self._traverse_loop(self.root, None)

        result = []
        for index in range(1, self.nbase + 1):
            x = self.BACKBONE_DISTANCE * float(self.bases[index].x)
            y = self.BACKBONE_DISTANCE * float(self.bases[index].y)
            if not math.isfinite(x) or not math.isfinite(y):
                raise RuntimeError("NAView generated a non-finite coordinate")
            result.append((x, y))
        return result

    def _read_in_bases(self, pair_table):
        self.bases[0].mate = 0
        self.bases[0].extracted = False
        pair_count = 0
        for index in range(1, self.nbase + 1):
            base = self.bases[index]
            base.extracted = False
            base.x = self.ANUM
            base.y = self.ANUM
            base.mate = int(pair_table[index])
            if base.mate > index:
                pair_count += 1
        if pair_count == 0:
            raise ValueError("NAView requires at least one planar base pair")

    def _find_regions(self):
        mark = [False] * (self.nbase + 1)
        self.nregion = 0
        index = 0
        while index <= self.nbase:
            mate = int(self.bases[index].mate)
            if mate != 0 and not mark[index]:
                region = self.regions[self.nregion]
                region.start1 = index
                region.end2 = mate
                mark[index] = True
                mark[mate] = True
                self.bases[index].region = region
                self.bases[mate].region = region

                index += 1
                mate -= 1
                while index < mate and int(self.bases[index].mate) == mate:
                    mark[index] = True
                    mark[mate] = True
                    self.bases[index].region = region
                    self.bases[mate].region = region
                    index += 1
                    mate -= 1

                index -= 1
                region.end1 = index
                region.start2 = mate + 1
                self.nregion += 1
            index += 1

    def _construct_loop(self, ibase):
        if self.loop_count >= len(self.loops):
            raise RuntimeError("NAView loop table overflow")
        result = self.loops[self.loop_count]
        self.loop_count += 1
        result.connections = []
        result.depth = 0
        result.radius = 0.0

        index = int(ibase)
        guard = 0
        while True:
            guard += 1
            if guard > (self.nbase + 2) * 4:
                raise RuntimeError("NAView loop construction did not terminate")

            mate = int(self.bases[index].mate)
            if mate != 0:
                region = self.bases[index].region
                if region is None:
                    raise RuntimeError("NAView base has no region")
                if not self.bases[region.start1].extracted:
                    if index == region.start1:
                        for item in (
                            region.start1,
                            region.end1,
                            region.start2,
                            region.end2,
                        ):
                            self.bases[item].extracted = True
                        child_loop = self._construct_loop(
                            region.end1 + 1 if region.end1 < self.nbase else 0
                        )
                    elif index == region.start2:
                        for item in (
                            region.start2,
                            region.end2,
                            region.start1,
                            region.end1,
                        ):
                            self.bases[item].extracted = True
                        child_loop = self._construct_loop(
                            region.end2 + 1 if region.end2 < self.nbase else 0
                        )
                    else:
                        raise RuntimeError("NAView loop construction invariant failed")

                    connection = DssrNaviewConnection()
                    connection.loop = child_loop
                    connection.region = region
                    if index == region.start1:
                        connection.start = region.start1
                        connection.end = region.end2
                    else:
                        connection.start = region.start2
                        connection.end = region.end1
                    result.connections.append(connection)

                    reverse = DssrNaviewConnection()
                    reverse.loop = result
                    reverse.region = region
                    if index == region.start1:
                        reverse.start = region.start2
                        reverse.end = region.end1
                    else:
                        reverse.start = region.start1
                        reverse.end = region.end2
                    child_loop.connections.append(reverse)

                index = mate

            index += 1
            if index > self.nbase:
                index = 0
            if index == ibase:
                break
        return result

    def _depth(self, loop):
        if loop.nconnection <= 1:
            return 0
        if loop.mark:
            return -1
        loop.mark = True
        count = 0
        result = 0
        for connection in loop.connections:
            depth = self._depth(connection.loop)
            if depth >= 0:
                count += 1
                if count == 1:
                    result = depth
                elif result > depth:
                    result = depth
        loop.mark = False
        return result + 1

    def _find_central_loop(self):
        active_loops = self.loops[: self.loop_count]
        for loop in active_loops:
            for candidate in active_loops:
                candidate.mark = False
            loop.depth = self._depth(loop)

        max_connections = 0
        max_depth = -1
        self.root = None
        for loop in active_loops:
            if loop.nconnection > max_connections or (
                loop.nconnection == max_connections and loop.depth > max_depth
            ):
                max_connections = loop.nconnection
                max_depth = loop.depth
                self.root = loop

    @staticmethod
    def _connected(connection, following):
        return bool(connection.extruded) or (
            int(connection.end) + 1 == int(following.start)
        )

    def _find_middle_connection(
        self, start, end, anchor_connection, anchor_in_loop, loop
    ):
        count = 0
        result = -1
        index = int(start)
        while True:
            count += 1
            if (
                anchor_connection is not None
                and loop.connections[index] is anchor_in_loop
            ):
                result = index
            if index == end:
                break
            index = (index + 1) % loop.nconnection
            if count > loop.nconnection * 3:
                raise RuntimeError("NAView connection walk did not terminate")

        if result == -1:
            index = int(start)
            for _unused in range(1, (count + 1) // 2):
                index = (index + 1) % loop.nconnection
            result = index
        return result

    def _determine_radius(self, loop, length_cutoff):
        minimum_radius = 0.7071068
        guard = 0
        while True:
            minimum_increment = 1.0e10
            numerator = 0.0
            denominator = 0.0
            minimum_index = 0

            for index, connection in enumerate(loop.connections):
                following = loop.connections[(index + 1) % loop.nconnection]
                end = int(connection.end)
                start = int(following.start)
                if start < end:
                    start += self.nbase + 1

                delta_angle = float(following.angle) - float(connection.angle)
                if delta_angle <= 0.0:
                    delta_angle += 2.0 * math.pi

                if not connection.extruded:
                    count = float(start - end)
                else:
                    count = 2.0 if delta_angle <= math.pi / 2.0 else 1.5
                if count <= 0.0:
                    count = 1.0e-6

                numerator += delta_angle * (1.0 / count + 1.0)
                denominator += delta_angle * delta_angle / count
                increment = delta_angle / count
                if (
                    increment < minimum_increment
                    and not connection.extruded
                    and count > 1.0
                ):
                    minimum_increment = increment
                    minimum_index = index

            radius = (
                numerator / denominator if denominator > 1.0e-14 else minimum_radius
            )
            radius = max(radius, minimum_radius)
            if minimum_increment * radius < length_cutoff:
                loop.connections[minimum_index].extruded = True
            else:
                break

            guard += 1
            if guard > loop.nconnection + 5:
                break

        if loop.radius > 0.0:
            radius = loop.radius
        else:
            loop.radius = radius

    def _traverse_loop(self, loop, anchor_connection):
        if loop.nconnection <= 0:
            return

        angle_increment = 2.0 * math.pi / float(self.nbase + 1)
        anchor_in_loop = None
        root_connection_index = -1

        for index, connection in enumerate(loop.connections):
            start_x = -math.sin(angle_increment * connection.start)
            start_y = math.cos(angle_increment * connection.start)
            end_x = -math.sin(angle_increment * connection.end)
            end_y = math.cos(angle_increment * connection.end)
            normal_x = end_y - start_y
            normal_y = start_x - end_x
            length = math.hypot(normal_x, normal_y)
            if length <= 1.0e-14:
                normal_x, normal_y, length = 1.0, 0.0, 1.0
            connection.xrad = normal_x / length
            connection.yrad = normal_y / length
            connection.angle = math.atan2(normal_y, normal_x) % (2.0 * math.pi)
            if (
                anchor_connection is not None
                and anchor_connection.region is connection.region
            ):
                anchor_in_loop = connection
                root_connection_index = index

        restart = True
        restart_count = 0
        while restart:
            restart = False
            restart_count += 1
            if restart_count > 100:
                raise RuntimeError("NAView exceeded its loop restart limit")

            self._determine_radius(loop, self.lencut)
            radius = loop.radius / self.RADIUS_REDUCTION_FACTOR
            if anchor_connection is None:
                center_x = 0.0
                center_y = 0.0
            else:
                origin_x = (
                    self.bases[anchor_in_loop.start].x
                    + self.bases[anchor_in_loop.end].x
                ) / 2.0
                origin_y = (
                    self.bases[anchor_in_loop.start].y
                    + self.bases[anchor_in_loop.end].y
                ) / 2.0
                center_x = origin_x - radius * anchor_in_loop.xrad
                center_y = origin_y - radius * anchor_in_loop.yrad

            connection_start = (
                0 if root_connection_index == -1 else root_connection_index
            )
            connection = loop.connections[connection_start]
            count = 0
            while True:
                previous_index = (connection_start - 1) % loop.nconnection
                previous = loop.connections[previous_index]
                if not self._connected(previous, connection):
                    break
                connection_start = previous_index
                connection = previous
                count += 1
                if count > loop.nconnection:
                    largest_angle = -1.0
                    largest_index = 0
                    for index, candidate in enumerate(loop.connections):
                        following = loop.connections[(index + 1) % loop.nconnection]
                        separation = (following.angle - candidate.angle) % (
                            2.0 * math.pi
                        )
                        if separation > largest_angle:
                            largest_angle = separation
                            largest_index = index
                    connection_end = largest_index
                    connection_start = (largest_index + 1) % loop.nconnection
                    break

            first_start = connection_start
            all_connections_done = False
            while not all_connections_done:
                count = 0
                connection_end = connection_start
                rooted = False
                while True:
                    connection = loop.connections[connection_end]
                    if connection_end == root_connection_index:
                        rooted = True
                    following_index = (connection_end + 1) % loop.nconnection
                    following = loop.connections[following_index]
                    if self._connected(connection, following):
                        count += 1
                        if count >= loop.nconnection:
                            break
                        connection_end = following_index
                    else:
                        break

                middle = self._find_middle_connection(
                    connection_start,
                    connection_end,
                    anchor_connection,
                    anchor_in_loop,
                    loop,
                )
                upward = middle
                downward = middle
                direction = 0

                while True:
                    if direction < 0:
                        current_index = upward
                    elif direction == 0:
                        current_index = middle
                    else:
                        current_index = downward

                    if current_index >= 0:
                        connection = loop.connections[current_index]
                        if (
                            anchor_connection is None
                            or anchor_in_loop is not connection
                        ):
                            if direction == 0:
                                half_angle = math.asin(min(1.0, 1.0 / (2.0 * radius)))
                                start_angle = connection.angle - half_angle
                                end_angle = connection.angle + half_angle
                                self.bases[connection.start].x = (
                                    center_x + radius * math.cos(start_angle)
                                )
                                self.bases[connection.start].y = (
                                    center_y + radius * math.sin(start_angle)
                                )
                                self.bases[connection.end].x = (
                                    center_x + radius * math.cos(end_angle)
                                )
                                self.bases[connection.end].y = (
                                    center_y + radius * math.sin(end_angle)
                                )
                            elif direction < 0:
                                following_index = (current_index + 1) % loop.nconnection
                                connection = loop.connections[current_index]
                                following = loop.connections[following_index]
                                angle = (connection.angle + following.angle) / 2.0
                                if connection.angle > following.angle:
                                    angle -= math.pi
                                line_x = math.sin(angle)
                                line_y = -math.cos(angle)
                                separation = (following.angle - connection.angle) % (
                                    2.0 * math.pi
                                )
                                multiplier = (
                                    2.0
                                    if connection.extruded
                                    and separation <= math.pi / 2.0
                                    else (1.5 if connection.extruded else 1.0)
                                )
                                self.bases[connection.end].x = (
                                    self.bases[following.start].x + multiplier * line_x
                                )
                                self.bases[connection.end].y = (
                                    self.bases[following.start].y + multiplier * line_y
                                )
                                self.bases[connection.start].x = (
                                    self.bases[connection.end].x + connection.yrad
                                )
                                self.bases[connection.start].y = (
                                    self.bases[connection.end].y - connection.xrad
                                )
                            else:
                                previous_index = (current_index - 1) % loop.nconnection
                                previous = loop.connections[previous_index]
                                connection = loop.connections[current_index]
                                angle = (previous.angle + connection.angle) / 2.0
                                if previous.angle > connection.angle:
                                    angle -= math.pi
                                line_x = -math.sin(angle)
                                line_y = math.cos(angle)
                                separation = (connection.angle - previous.angle) % (
                                    2.0 * math.pi
                                )
                                multiplier = (
                                    2.0
                                    if previous.extruded and separation <= math.pi / 2.0
                                    else (1.5 if previous.extruded else 1.0)
                                )
                                self.bases[connection.start].x = (
                                    self.bases[previous.end].x + multiplier * line_x
                                )
                                self.bases[connection.start].y = (
                                    self.bases[previous.end].y + multiplier * line_y
                                )
                                self.bases[connection.end].x = (
                                    self.bases[connection.start].x - connection.yrad
                                )
                                self.bases[connection.end].y = (
                                    self.bases[connection.start].y + connection.xrad
                                )

                    if direction < 0:
                        if downward == connection_end:
                            downward = -1
                        elif downward >= 0:
                            downward = (downward + 1) % loop.nconnection
                        direction = 1
                    else:
                        if upward == connection_start:
                            upward = -1
                        elif upward >= 0:
                            upward = (upward - 1) % loop.nconnection
                        direction = -1
                    if upward == -1 and downward == -1:
                        break

                next_start = (connection_end + 1) % loop.nconnection
                if connection_end != connection_start and not (
                    connection_start == first_start and next_start == first_start
                ):
                    first_connection = loop.connections[connection_start]
                    last_connection = loop.connections[connection_end]
                    delta_x = (
                        self.bases[last_connection.end].x
                        - self.bases[first_connection.start].x
                    )
                    delta_y = (
                        self.bases[last_connection.end].y
                        - self.bases[first_connection.start].y
                    )
                    middle_x = self.bases[first_connection.start].x + delta_x / 2.0
                    middle_y = self.bases[first_connection.start].y + delta_y / 2.0
                    length = math.hypot(delta_x, delta_y)
                    if length > 1.0e-12:
                        tangent_x = delta_x / length
                        tangent_y = delta_y / length
                        vector_x = (center_x - middle_x) / length
                        vector_y = (center_y - middle_y) / length
                        dot = vector_x * tangent_x + vector_y * tangent_y
                        normal_x = dot * tangent_x - vector_x
                        normal_y = dot * tangent_y - vector_y
                        normal_length = math.hypot(normal_x, normal_y)
                        if normal_length > 1.0e-12:
                            normal_x /= normal_length
                            normal_y /= normal_length
                            start_angle = math.atan2(
                                self.bases[first_connection.start].y - center_y,
                                self.bases[first_connection.start].x - center_x,
                            ) % (2.0 * math.pi)
                            end_angle = math.atan2(
                                self.bases[last_connection.end].y - center_y,
                                self.bases[last_connection.end].x - center_x,
                            ) % (2.0 * math.pi)
                            if end_angle < start_angle:
                                end_angle += 2.0 * math.pi
                            sign = -1.0 if end_angle - start_angle > math.pi else 1.0
                            new_middle_x = center_x + sign * radius * normal_x
                            new_middle_y = center_y + sign * radius * normal_y
                            if rooted:
                                center_x -= new_middle_x - middle_x
                                center_y -= new_middle_y - middle_y
                            else:
                                current_index = connection_start
                                while True:
                                    candidate = loop.connections[current_index]
                                    for base_index in (
                                        candidate.start,
                                        candidate.end,
                                    ):
                                        self.bases[base_index].x += (
                                            new_middle_x - middle_x
                                        )
                                        self.bases[base_index].y += (
                                            new_middle_y - middle_y
                                        )
                                    if current_index == connection_end:
                                        break
                                    current_index = (
                                        current_index + 1
                                    ) % loop.nconnection

                connection_start = next_start
                all_connections_done = connection_start == first_start

            for index, connection in enumerate(loop.connections):
                following = loop.connections[(index + 1) % loop.nconnection]
                delta_x = self.bases[connection.end].x - center_x
                delta_y = self.bases[connection.end].y - center_y
                radius_current = math.hypot(delta_x, delta_y)
                angle_current = math.atan2(delta_y, delta_x) % (2.0 * math.pi)

                delta_x = self.bases[following.start].x - center_x
                delta_y = self.bases[following.start].y - center_y
                radius_following = math.hypot(delta_x, delta_y)
                angle_following = math.atan2(delta_y, delta_x) % (2.0 * math.pi)
                if angle_following < angle_current:
                    angle_following += 2.0 * math.pi

                sweep = angle_following - angle_current
                expected = (following.angle - connection.angle) % (2.0 * math.pi)
                if abs(sweep - expected) > math.pi:
                    if (
                        not connection.extruded
                        and (following.start - connection.end) != 1
                    ):
                        connection.extruded = True
                        restart = True
                        break

                if connection.extruded:
                    self._construct_extruded_segment(connection, following)
                else:
                    count = int(following.start) - int(connection.end)
                    if count < 0:
                        count += self.nbase + 1
                    if count > 0:
                        increment = sweep / float(count)
                        for offset in range(1, count):
                            base_index = int(connection.end) + offset
                            if base_index > self.nbase:
                                base_index -= self.nbase + 1
                            angle = angle_current + offset * increment
                            if abs(sweep) > 1.0e-12:
                                local_radius = (
                                    radius_current
                                    + (radius_following - radius_current)
                                    * (angle - angle_current)
                                    / sweep
                                )
                            else:
                                local_radius = radius_current
                            self.bases[base_index].x = (
                                center_x + local_radius * math.cos(angle)
                            )
                            self.bases[base_index].y = (
                                center_y + local_radius * math.sin(angle)
                            )

            if restart:
                continue

        for index, connection in enumerate(loop.connections):
            if root_connection_index != index:
                self._generate_region(connection)
                self._traverse_loop(connection.loop, connection)

        count = 0
        sum_x = 0.0
        sum_y = 0.0
        for index, connection in enumerate(loop.connections):
            following = loop.connections[(index + 1) % loop.nconnection]
            count += 2
            sum_x += self.bases[connection.start].x + self.bases[connection.end].x
            sum_y += self.bases[connection.start].y + self.bases[connection.end].y
            if not connection.extruded:
                base_index = int(connection.end) + 1
                while base_index != int(following.start):
                    if base_index > self.nbase:
                        base_index -= self.nbase + 1
                    count += 1
                    sum_x += self.bases[base_index].x
                    sum_y += self.bases[base_index].y
                    base_index += 1

        if count > 0:
            loop.x = sum_x / float(count)
            loop.y = sum_y / float(count)

    def _generate_region(self, connection):
        region = connection.region
        if connection.start == region.start1:
            start = region.start1
            end = region.end1
        else:
            start = region.start2
            end = region.end2

        length = 0
        for base_index in range(start + 1, end + 1):
            length += 1
            self.bases[base_index].x = (
                self.bases[connection.start].x
                + self.HELIX_FACTOR * length * connection.xrad
            )
            self.bases[base_index].y = (
                self.bases[connection.start].y
                + self.HELIX_FACTOR * length * connection.yrad
            )
            mate = int(self.bases[base_index].mate)
            self.bases[mate].x = (
                self.bases[connection.end].x
                + self.HELIX_FACTOR * length * connection.xrad
            )
            self.bases[mate].y = (
                self.bases[connection.end].y
                + self.HELIX_FACTOR * length * connection.yrad
            )

    def _construct_extruded_segment(self, connection, following):
        start_angle = float(connection.angle)
        end_angle_1 = float(following.angle)
        end_angle_2 = end_angle_1
        if end_angle_2 < start_angle:
            end_angle_2 += 2.0 * math.pi
        average_angle = (start_angle + end_angle_2) / 2.0

        start = int(connection.end)
        end = int(following.start)
        count = end - start
        if count < 0:
            count += self.nbase + 1
        separation = (following.angle - connection.angle) % (2.0 * math.pi)

        if count == 2:
            self._construct_circle_segment(start, end)
            return

        delta_x = self.bases[end].x - self.bases[start].x
        delta_y = self.bases[end].y - self.bases[start].y
        distance = math.hypot(delta_x, delta_y)
        if distance <= 1.0e-12:
            return
        delta_x /= distance
        delta_y /= distance

        if distance >= 1.5 and separation <= math.pi / 2.0:
            next_start = start + 1
            if next_start > self.nbase:
                next_start -= self.nbase + 1
            previous_end = end - 1
            if previous_end < 0:
                previous_end += self.nbase + 1
            self.bases[next_start].x = self.bases[start].x + 0.5 * delta_x
            self.bases[next_start].y = self.bases[start].y + 0.5 * delta_y
            self.bases[previous_end].x = self.bases[end].x - 0.5 * delta_x
            self.bases[previous_end].y = self.bases[end].y - 0.5 * delta_y
            start = next_start
            end = previous_end

        collision = True
        while collision and count > 1:
            collision = False
            self._construct_circle_segment(start, end)

            next_start = start + 1
            if next_start > self.nbase:
                next_start -= self.nbase + 1
            angle_1 = math.atan2(
                self.bases[next_start].y - self.bases[start].y,
                self.bases[next_start].x - self.bases[start].x,
            ) % (2.0 * math.pi)
            if (angle_1 - start_angle) % (2.0 * math.pi) > math.pi:
                collision = True

            previous_end = end - 1
            if previous_end < 0:
                previous_end += self.nbase + 1
            angle_2 = math.atan2(
                self.bases[previous_end].y - self.bases[end].y,
                self.bases[previous_end].x - self.bases[end].x,
            ) % (2.0 * math.pi)
            if (end_angle_1 - angle_2) % (2.0 * math.pi) > math.pi:
                collision = True

            if collision:
                angle = min(average_angle, start_angle + 0.5)
                self.bases[next_start].x = self.bases[start].x + math.cos(angle)
                self.bases[next_start].y = self.bases[start].y + math.sin(angle)
                start = next_start

                angle = max(average_angle, end_angle_2 - 0.5)
                self.bases[previous_end].x = self.bases[end].x + math.cos(angle)
                self.bases[previous_end].y = self.bases[end].y + math.sin(angle)
                end = previous_end
                count -= 2

    def _construct_circle_segment(self, start, end):
        delta_x = self.bases[end].x - self.bases[start].x
        delta_y = self.bases[end].y - self.bases[start].y
        distance = math.hypot(delta_x, delta_y)
        length = end - start
        if length < 0:
            length += self.nbase + 1
        if length <= 0:
            return

        if distance >= length:
            if distance <= 1.0e-12:
                return
            delta_x /= distance
            delta_y /= distance
            for offset in range(1, length):
                base_index = start + offset
                if base_index > self.nbase:
                    base_index -= self.nbase + 1
                self.bases[base_index].x = self.bases[
                    start
                ].x + delta_x * offset / float(length)
                self.bases[base_index].y = self.bases[
                    start
                ].y + delta_y * offset / float(length)
            return

        self._find_center_for_arc(length - 1, distance)
        if distance <= 1.0e-12:
            return
        delta_x /= distance
        delta_y /= distance
        middle_x = self.bases[start].x + delta_x * distance / 2.0
        middle_y = self.bases[start].y + delta_y * distance / 2.0
        normal_x = delta_y
        normal_y = -delta_x
        center_x = middle_x + self._h * normal_x
        center_y = middle_y + self._h * normal_y
        vector_x = self.bases[start].x - center_x
        vector_y = self.bases[start].y - center_y
        radius = math.hypot(vector_x, vector_y)
        angle = math.atan2(vector_y, vector_x)

        for offset in range(1, length):
            base_index = start + offset
            if base_index > self.nbase:
                base_index -= self.nbase + 1
            self.bases[base_index].x = center_x + radius * math.cos(
                angle + offset * self.angleinc
            )
            self.bases[base_index].y = center_y + radius * math.sin(
                angle + offset * self.angleinc
            )

    def _find_center_for_arc(self, count, chord):
        upper = (count + 1.0) / math.pi
        lower = -upper - chord / (count + 1.000001 - chord)
        if chord < 1.0:
            lower = 0.0

        height = 0.0
        theta = 0.0
        error = 0.0
        for _iteration in range(self.MAXITER):
            height = (upper + lower) / 2.0
            radius = math.sqrt(height * height + chord * chord / 4.0)
            if radius <= 1.0e-14:
                break
            discriminant = 1.0 - 0.5 / (radius * radius)
            discriminant = max(-1.0, min(1.0, discriminant))
            theta = math.acos(discriminant)
            phi = math.acos(max(-1.0, min(1.0, height / radius)))
            error = theta * (count + 1) + 2.0 * phi - 2.0 * math.pi
            if error > 0.0:
                lower = height
            else:
                upper = height
            if abs(error) <= 0.0001:
                break

        self._h = height
        self.angleinc = theta


# RNA layout algorithms
class Dssr2DLayout:
    """Dependency-free layout algorithms for RNA graphs."""

    @staticmethod
    def circular(model):
        n = len(model.nts)
        if n <= 0:
            return []
        break_gap = 0.40
        total_gap = break_gap * len(model.chain_breaks)
        step = (2.0 * math.pi - total_gap) / max(1, n)
        radius = max(150.0, (n * 38.0) / (2.0 * math.pi))
        theta = -math.pi / 2.0
        out = []
        for i in range(n):
            out.append((radius * math.cos(theta), radius * math.sin(theta)))
            theta += step
            if i in model.chain_breaks:
                theta += break_gap
        return out

    @staticmethod
    def linear(model):
        out = []
        x = 0.0
        for i in range(len(model.nts)):
            out.append((x, 0.0))
            x += 42.0
            if i in model.chain_breaks:
                x += 70.0
        if out:
            mid = (out[0][0] + out[-1][0]) / 2.0
            out = [(x - mid, y) for x, y in out]
        return out

    @staticmethod
    def _v_add(a, b):
        return (a[0] + b[0], a[1] + b[1])

    @staticmethod
    def _v_mul(a, scalar):
        return (a[0] * scalar, a[1] * scalar)

    @staticmethod
    def _v_norm(a):
        length = math.hypot(a[0], a[1])
        if length <= 1.0e-12:
            return (0.0, 1.0)
        return (a[0] / length, a[1] / length)

    @staticmethod
    def _v_rotate(a, angle):
        c = math.cos(angle)
        s = math.sin(angle)
        return (a[0] * c - a[1] * s, a[0] * s + a[1] * c)

    @staticmethod
    def _planar_pair_table(model):
        n = len(model.nts)
        table = [-1] * n
        for pair in model.planar_secondary_pairs():
            i = int(pair.get("i", -1))
            j = int(pair.get("j", -1))
            if i > j:
                i, j = j, i
            if i < 0 or j >= n or i == j:
                continue
            if table[i] == -1 and table[j] == -1:
                table[i] = j
                table[j] = i
        return table

    @staticmethod
    def _stem_tree(pair_table, chain_breaks):
        stems = []
        n = len(pair_table)
        breaks = set(chain_breaks)

        for i in range(n):
            j = pair_table[i]
            if j <= i:
                continue

            previous_is_same_stem = (
                i > 0
                and (i - 1) not in breaks
                and j < n - 1
                and j not in breaks
                and pair_table[i - 1] == j + 1
            )
            if previous_is_same_stem:
                continue

            pairs = []
            a, b = i, j
            while a < b and pair_table[a] == b:
                pairs.append((a, b))
                if a in breaks or (b - 1) in breaks:
                    break
                a += 1
                b -= 1

            if not pairs:
                continue
            stems.append(
                {
                    "outer_i": pairs[0][0],
                    "outer_j": pairs[0][1],
                    "inner_i": pairs[-1][0],
                    "inner_j": pairs[-1][1],
                    "pairs": pairs,
                    "parent": None,
                    "children": [],
                }
            )

        # The smallest containing stem is the direct parent.
        for stem in stems:
            containers = [
                other
                for other in stems
                if other is not stem
                and other["inner_i"] < stem["outer_i"]
                and stem["outer_j"] < other["inner_j"]
            ]
            if containers:
                parent = min(
                    containers,
                    key=lambda other: other["inner_j"] - other["inner_i"],
                )
                stem["parent"] = parent
                parent["children"].append(stem)

        for stem in stems:
            stem["children"].sort(key=lambda child: child["outer_i"])
        return stems

    @staticmethod
    def _radiate_general(model):
        """Draw a planar stem/loop scaffold and overlay pseudoknots later."""
        n = len(model.nts)
        if n <= 0:
            return []

        pair_table = Dssr2DLayout._planar_pair_table(model)
        stems = Dssr2DLayout._stem_tree(pair_table, model.chain_breaks)
        if not stems:
            return Dssr2DLayout.circular(model)

        positions = [None] * n
        directions = [None] * n
        rise = 40.0
        half_width = 18.0
        placed_stems = set()

        def place_stem(stem, outer_center, direction):
            marker = id(stem)
            if marker in placed_stems:
                return
            placed_stems.add(marker)

            direction = Dssr2DLayout._v_norm(direction)
            normal = (-direction[1], direction[0])
            pairs = stem["pairs"]

            for k, (i, j) in enumerate(pairs):
                center = Dssr2DLayout._v_add(
                    outer_center,
                    Dssr2DLayout._v_mul(direction, k * rise),
                )
                positions[i] = Dssr2DLayout._v_add(
                    center,
                    Dssr2DLayout._v_mul(normal, half_width),
                )
                positions[j] = Dssr2DLayout._v_add(
                    center,
                    Dssr2DLayout._v_mul(normal, -half_width),
                )
                directions[i] = direction
                directions[j] = direction

            inner_center = Dssr2DLayout._v_add(
                outer_center,
                Dssr2DLayout._v_mul(direction, (len(pairs) - 1) * rise),
            )
            children = list(stem["children"])
            count = len(children)
            if count <= 0:
                return

            if count == 1:
                angles = [0.0]
            else:
                max_angle = math.radians(min(78.0, 35.0 + 15.0 * (count - 1)))
                angles = [
                    -max_angle + (2.0 * max_angle * k / float(count - 1))
                    for k in range(count)
                ]

            branch_distance = 85.0 + 10.0 * max(0, count - 2)
            for child, angle in zip(children, angles):
                child_direction = Dssr2DLayout._v_rotate(direction, angle)
                child_center = Dssr2DLayout._v_add(
                    inner_center,
                    Dssr2DLayout._v_mul(child_direction, branch_distance),
                )
                place_stem(child, child_center, child_direction)

        top_level = [stem for stem in stems if stem["parent"] is None]
        top_level.sort(key=lambda stem: stem["outer_i"])
        top_count = len(top_level)

        if top_count == 1:
            place_stem(top_level[0], (0.0, 0.0), (0.0, 1.0))
        elif top_count <= 4:
            max_angle = math.radians(70.0)
            angles = [
                -max_angle + (2.0 * max_angle * k / float(top_count - 1))
                for k in range(top_count)
            ]
            for stem, angle in zip(top_level, angles):
                direction = Dssr2DLayout._v_rotate((0.0, 1.0), angle)
                center = Dssr2DLayout._v_mul(direction, 80.0)
                place_stem(stem, center, direction)
        else:
            for k, stem in enumerate(top_level):
                angle = -math.pi / 2.0 + 2.0 * math.pi * k / float(top_count)
                direction = (math.cos(angle), math.sin(angle))
                center = Dssr2DLayout._v_mul(direction, 120.0)
                place_stem(stem, center, direction)

        def fill_unknown_run(start, end, segment_start, segment_end):
            count = end - start + 1
            previous = start - 1 if start > segment_start else None
            following = end + 1 if end < segment_end else None

            if previous is not None and following is not None:
                a = positions[previous]
                b = positions[following]
                if a is None or b is None:
                    return

                average_direction = (0.0, 0.0)
                for value in (directions[previous], directions[following]):
                    if value is not None:
                        average_direction = Dssr2DLayout._v_add(
                            average_direction, value
                        )
                if math.hypot(average_direction[0], average_direction[1]) < 0.1:
                    chord = (b[0] - a[0], b[1] - a[1])
                    average_direction = (-chord[1], chord[0])
                average_direction = Dssr2DLayout._v_norm(average_direction)

                if pair_table[previous] == following:
                    midpoint = ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)
                    chord_length = math.hypot(b[0] - a[0], b[1] - a[1])
                    half_chord = chord_length * 0.5
                    radius = max(
                        half_chord + 1.0,
                        ((count + 1) * 34.0 + chord_length) / (2.0 * math.pi),
                    )
                    center_distance = math.sqrt(
                        max(1.0, radius * radius - half_chord * half_chord)
                    )
                    circle_center = Dssr2DLayout._v_add(
                        midpoint,
                        Dssr2DLayout._v_mul(average_direction, center_distance),
                    )
                    theta_a = math.atan2(
                        a[1] - circle_center[1], a[0] - circle_center[0]
                    )
                    theta_b = math.atan2(
                        b[1] - circle_center[1], b[0] - circle_center[0]
                    )
                    minor_sweep = (theta_b - theta_a) % (2.0 * math.pi)
                    long_sweep = 2.0 * math.pi - minor_sweep
                    for offset, index in enumerate(range(start, end + 1), 1):
                        t = offset / float(count + 1)
                        theta = theta_a - long_sweep * t
                        positions[index] = (
                            circle_center[0] + radius * math.cos(theta),
                            circle_center[1] + radius * math.sin(theta),
                        )
                        directions[index] = average_direction
                    return

                height = max(35.0, min(150.0, 18.0 * count + 10.0))
                midpoint = ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)
                control = Dssr2DLayout._v_add(
                    midpoint,
                    Dssr2DLayout._v_mul(average_direction, height),
                )
                for offset, index in enumerate(range(start, end + 1), 1):
                    t = offset / float(count + 1)
                    one_minus = 1.0 - t
                    positions[index] = (
                        one_minus * one_minus * a[0]
                        + 2.0 * one_minus * t * control[0]
                        + t * t * b[0],
                        one_minus * one_minus * a[1]
                        + 2.0 * one_minus * t * control[1]
                        + t * t * b[1],
                    )
                    directions[index] = average_direction
                return

            if following is not None and positions[following] is not None:
                anchor = positions[following]
                direction = directions[following] or (0.0, 1.0)
                normal = (-direction[1], direction[0])
                for offset, index in enumerate(range(end, start - 1, -1), 1):
                    positions[index] = Dssr2DLayout._v_add(
                        anchor,
                        Dssr2DLayout._v_add(
                            Dssr2DLayout._v_mul(direction, -36.0 * offset),
                            Dssr2DLayout._v_mul(normal, -10.0 * offset),
                        ),
                    )
                    directions[index] = direction
                return

            if previous is not None and positions[previous] is not None:
                anchor = positions[previous]
                direction = directions[previous] or (0.0, 1.0)
                normal = (-direction[1], direction[0])
                for offset, index in enumerate(range(start, end + 1), 1):
                    positions[index] = Dssr2DLayout._v_add(
                        anchor,
                        Dssr2DLayout._v_add(
                            Dssr2DLayout._v_mul(direction, -36.0 * offset),
                            Dssr2DLayout._v_mul(normal, 10.0 * offset),
                        ),
                    )
                    directions[index] = direction
                return

            for offset, index in enumerate(range(start, end + 1)):
                positions[index] = (offset * 40.0, 0.0)
                directions[index] = (1.0, 0.0)

        segments = []
        segment_start = 0
        for break_after in sorted(model.chain_breaks):
            if break_after >= segment_start:
                segments.append((segment_start, min(n - 1, break_after)))
                segment_start = break_after + 1
        if segment_start < n:
            segments.append((segment_start, n - 1))

        for segment_start, segment_end in segments:
            index = segment_start
            while index <= segment_end:
                if positions[index] is not None:
                    index += 1
                    continue
                run_start = index
                while index <= segment_end and positions[index] is None:
                    index += 1
                fill_unknown_run(
                    run_start,
                    index - 1,
                    segment_start,
                    segment_end,
                )

        fallback = Dssr2DLayout.circular(model)
        for i in range(n):
            if positions[i] is None:
                positions[i] = fallback[i]

        center_x = sum(point[0] for point in positions) / float(n)
        center_y = sum(point[1] for point in positions) / float(n)
        return [(point[0] - center_x, point[1] - center_y) for point in positions]

    @staticmethod
    def _solve_circle(edge_lengths):
        """Solve a circle whose successive chord lengths close one revolution."""
        lengths = [max(1.0, float(value)) for value in edge_lengths]
        if not lengths:
            return 30.0, []

        lower = max(lengths) * 0.5 + 1.0e-7

        def angles_at(radius):
            return [
                2.0 * math.asin(min(1.0, length / (2.0 * radius))) for length in lengths
            ]

        lower_angles = angles_at(lower)
        lower_total = sum(lower_angles)
        if lower_total < 2.0 * math.pi:
            scale = (2.0 * math.pi) / max(1.0e-12, lower_total)
            return lower, [angle * scale for angle in lower_angles]

        upper = max(sum(lengths), lower * 2.0)
        while sum(angles_at(upper)) > 2.0 * math.pi:
            upper *= 2.0

        for _ in range(80):
            middle = 0.5 * (lower + upper)
            if sum(angles_at(middle)) > 2.0 * math.pi:
                lower = middle
            else:
                upper = middle

        radius = 0.5 * (lower + upper)
        return radius, angles_at(radius)

    @staticmethod
    def _normalize(vector, fallback=(0.0, 1.0)):
        length = math.hypot(float(vector[0]), float(vector[1]))
        if length <= 1.0e-12:
            return fallback
        return (float(vector[0]) / length, float(vector[1]) / length)

    @staticmethod
    def _quadratic_equal(indices, start, stop, control, positions):
        """Place indices at near-equal arc-length intervals on a quadratic curve."""
        indices = list(indices)
        if not indices:
            return

        samples = max(160, 28 * (len(indices) + 1))
        points = []
        for sample in range(samples + 1):
            t = sample / float(samples)
            u = 1.0 - t
            points.append(
                (
                    u * u * start[0] + 2.0 * u * t * control[0] + t * t * stop[0],
                    u * u * start[1] + 2.0 * u * t * control[1] + t * t * stop[1],
                )
            )

        cumulative = [0.0]
        for first, second in zip(points, points[1:]):
            cumulative.append(
                cumulative[-1] + math.hypot(second[0] - first[0], second[1] - first[1])
            )
        total = cumulative[-1]
        if total <= 1.0e-12:
            return

        for offset, index in enumerate(indices, 1):
            target = total * offset / float(len(indices) + 1)
            right = bisect.bisect_left(cumulative, target)
            right = min(max(1, right), len(points) - 1)
            left = right - 1
            span = cumulative[right] - cumulative[left]
            fraction = 0.0 if span <= 1.0e-12 else (target - cumulative[left]) / span
            positions[index] = (
                points[left][0] + fraction * (points[right][0] - points[left][0]),
                points[left][1] + fraction * (points[right][1] - points[left][1]),
            )

    @staticmethod
    def _place_loop_circle(
        indices,
        start,
        stop,
        outward_direction,
        positions,
        backbone_distance=34.0,
        pair_distance=42.0,
    ):
        """Place a hairpin on the long arc opposite its supporting base pair."""
        indices = list(indices)
        if not indices:
            return

        outward = Dssr2DLayout._normalize(outward_direction)
        radius, angles = Dssr2DLayout._solve_circle(
            [float(backbone_distance)] * (len(indices) + 1) + [float(pair_distance)]
        )

        midpoint = (
            0.5 * (start[0] + stop[0]),
            0.5 * (start[1] + stop[1]),
        )
        half_chord = 0.5 * math.hypot(stop[0] - start[0], stop[1] - start[1])
        center_distance = math.sqrt(max(0.0, radius * radius - half_chord * half_chord))
        center = (
            midpoint[0] + outward[0] * center_distance,
            midpoint[1] + outward[1] * center_distance,
        )
        theta_start = math.atan2(start[1] - center[1], start[0] - center[0])

        candidates = []
        for sign in (-1.0, 1.0):
            theta = theta_start
            candidate = []
            for offset, index in enumerate(indices):
                theta += sign * angles[offset]
                candidate.append(
                    (
                        center[0] + radius * math.cos(theta),
                        center[1] + radius * math.sin(theta),
                    )
                )
            score = sum(
                (point[0] - midpoint[0]) * outward[0]
                + (point[1] - midpoint[1]) * outward[1]
                for point in candidate
            ) / float(max(1, len(candidate)))
            candidates.append((score, candidate))

        chosen = max(candidates, key=lambda item: item[0])[1]
        for index, point in zip(indices, chosen):
            positions[index] = point

    @staticmethod
    def _detect_trna_topology(model):
        """Detect a tRNA-like cloverleaf from topology, never from a PDB name."""
        total = len(model.nts)
        if total < 55 or total > 110 or model.chain_count() != 1:
            return None

        pair_table = Dssr2DLayout._planar_pair_table(model)
        stems = Dssr2DLayout._stem_tree(pair_table, model.chain_breaks)
        candidates = []
        for stem in stems:
            children = list(stem.get("children", []))
            if len(stem.get("pairs", [])) < 5 or len(children) != 3:
                continue
            if any(len(child.get("pairs", [])) < 3 for child in children):
                continue
            if stem.get("outer_i", 999999) > 5:
                continue
            if (total - 1) - stem.get("outer_j", -1) > 10:
                continue
            ordered = sorted(children, key=lambda child: child["outer_i"])
            if any(child["inner_j"] - child["inner_i"] < 4 for child in ordered):
                continue
            score = (
                20 * len(stem.get("pairs", []))
                + sum(len(child.get("pairs", [])) for child in ordered)
                - stem["outer_i"]
                - ((total - 1) - stem["outer_j"])
            )
            candidates.append((score, stem, ordered))

        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0], reverse=True)
        _score, root, arms = candidates[0]
        return {
            "root": root,
            "arms": arms,
            "pair_table": pair_table,
            "stems": stems,
        }

    @staticmethod
    def _trna_cloverleaf(model, topology):
        """Generate a clean, conventional four-arm tRNA cloverleaf."""
        total = len(model.nts)
        root = topology["root"]
        arms = topology["arms"]
        positions = [None] * total

        pair_distance = 42.0
        helix_rise = 38.0
        backbone_distance = 34.0

        def place_stem(stem, outer_center, direction):
            direction = Dssr2DLayout._normalize(direction)
            normal = (-direction[1], direction[0])
            for offset, (left, right) in enumerate(stem["pairs"]):
                center = (
                    outer_center[0] + direction[0] * helix_rise * offset,
                    outer_center[1] + direction[1] * helix_rise * offset,
                )
                positions[left] = (
                    center[0] + normal[0] * pair_distance * 0.5,
                    center[1] + normal[1] * pair_distance * 0.5,
                )
                positions[right] = (
                    center[0] - normal[0] * pair_distance * 0.5,
                    center[1] - normal[1] * pair_distance * 0.5,
                )

        root_inner_center = (0.0, -88.0)
        root_direction = (0.0, 1.0)
        root_outer_center = (
            root_inner_center[0],
            root_inner_center[1] - helix_rise * (len(root["pairs"]) - 1),
        )
        place_stem(root, root_outer_center, root_direction)

        arm_centers = [(-150.0, -4.0), (0.0, 160.0), (150.0, -4.0)]
        arm_directions = [(-1.0, 0.0), (0.0, 1.0), (1.0, 0.0)]

        for arm, center, direction in zip(arms, arm_centers, arm_directions):
            place_stem(arm, center, direction)
            loop_indices = range(arm["inner_i"] + 1, arm["inner_j"])
            Dssr2DLayout._place_loop_circle(
                loop_indices,
                positions[arm["inner_i"]],
                positions[arm["inner_j"]],
                direction,
                positions,
                backbone_distance=backbone_distance,
                pair_distance=pair_distance,
            )

        junctions = [
            (root["inner_i"], arms[0]["outer_i"], (-106.0, -106.0)),
            (arms[0]["outer_j"], arms[1]["outer_i"], (-122.0, 96.0)),
            (arms[1]["outer_j"], arms[2]["outer_i"], (122.0, 96.0)),
            (arms[2]["outer_j"], root["inner_j"], (106.0, -106.0)),
        ]
        for first, last, control in junctions:
            Dssr2DLayout._quadratic_equal(
                range(first + 1, last),
                positions[first],
                positions[last],
                control,
                positions,
            )

        for offset, index in enumerate(range(root["outer_i"] - 1, -1, -1), 1):
            anchor = positions[root["outer_i"]]
            positions[index] = (
                anchor[0] - backbone_distance * offset,
                anchor[1] - 8.0 * offset,
            )
        for offset, index in enumerate(range(root["outer_j"] + 1, total), 1):
            anchor = positions[root["outer_j"]]
            positions[index] = (
                anchor[0] + backbone_distance * offset,
                anchor[1] - 8.0 * offset,
            )

        fallback = Dssr2DLayout._radiate_general(model)
        for index in range(total):
            if positions[index] is None:
                positions[index] = fallback[index]

        center_x = sum(point[0] for point in positions) / float(total)
        center_y = sum(point[1] for point in positions) / float(total)
        model._dssr2d_layout_variant = "tRNA cloverleaf"
        return [(point[0] - center_x, point[1] - center_y) for point in positions]

    @staticmethod
    def _pair_table(model, start=0, end=None):
        """Create a 1-based NAView pair table from the planar DSSR scaffold."""
        total = len(model.nts)
        if end is None:
            end = total - 1
        start = max(0, int(start))
        end = min(total - 1, int(end))
        count = max(0, end - start + 1)
        table = [count] + [0] * count
        for pair in model.planar_secondary_pairs():
            first = int(pair.get("i", -1))
            second = int(pair.get("j", -1))
            if first > second:
                first, second = second, first
            if first < start or second > end:
                continue
            local_first = first - start + 1
            local_second = second - start + 1
            if table[local_first] == 0 and table[local_second] == 0:
                table[local_first] = local_second
                table[local_second] = local_first
        return table

    @staticmethod
    def _center(points):
        if not points:
            return []
        center_x = sum(point[0] for point in points) / float(len(points))
        center_y = sum(point[1] for point in points) / float(len(points))
        return [(point[0] - center_x, point[1] - center_y) for point in points]

    @staticmethod
    def _rotate(points, angle):
        cosine = math.cos(angle)
        sine = math.sin(angle)
        return [
            (
                point[0] * cosine - point[1] * sine,
                point[0] * sine + point[1] * cosine,
            )
            for point in points
        ]

    @staticmethod
    def _standardize_trna_orientation(model, points):
        """Orient tRNA with acceptor stem down, anticodon arm up."""
        try:
            topology = Dssr2DLayout._detect_trna_topology(model)
        except Exception:
            topology = None
        if not topology or len(points) != len(model.nts):
            return points, False

        root = topology["root"]
        arms = topology["arms"]
        outer_first, outer_second = root["pairs"][0]
        inner_first, inner_second = root["pairs"][-1]
        outer_center = (
            0.5 * (points[outer_first][0] + points[outer_second][0]),
            0.5 * (points[outer_first][1] + points[outer_second][1]),
        )
        inner_center = (
            0.5 * (points[inner_first][0] + points[inner_second][0]),
            0.5 * (points[inner_first][1] + points[inner_second][1]),
        )
        vector_x = outer_center[0] - inner_center[0]
        vector_y = outer_center[1] - inner_center[1]
        if math.hypot(vector_x, vector_y) > 1.0e-8:
            current_angle = math.atan2(vector_y, vector_x)
            points = Dssr2DLayout._rotate(points, math.pi / 2.0 - current_angle)

        if len(arms) >= 3:
            d_indices = list(
                range(int(arms[0]["outer_i"]), int(arms[0]["outer_j"]) + 1)
            )
            t_indices = list(
                range(int(arms[2]["outer_i"]), int(arms[2]["outer_j"]) + 1)
            )
            d_x = sum(points[index][0] for index in d_indices) / float(
                max(1, len(d_indices))
            )
            t_x = sum(points[index][0] for index in t_indices) / float(
                max(1, len(t_indices))
            )
            if d_x > t_x:
                points = [(-point[0], point[1]) for point in points]

            v_start, v_end = int(arms[1]["outer_j"]) + 1, int(arms[2]["outer_i"])
            if 3 <= v_end - v_start <= 7:
                p1, p2 = points[v_start - 1], points[v_end]
                ctrl = (max(p1[0], p2[0]) + 48.0, min(p1[1], p2[1]) - 68.0)
                Dssr2DLayout._quadratic_equal(
                    range(v_start, v_end), p1, p2, ctrl, points
                )

        return Dssr2DLayout._center(points), True

    @staticmethod
    def _naview_layout(model):
        """Return one coherent NAView-style scientific layout."""
        total = len(model.nts)
        if total <= 0:
            return []

        scaffold_pair_count = len(model.planar_secondary_pairs())
        if scaffold_pair_count <= 0:
            model._dssr2d_layout_variant = "unpaired circular"
            return Dssr2DLayout.circular(model)

        if model.chain_count() != 1:
            model._dssr2d_layout_variant = "multi-chain radiate fallback"
            return Dssr2DLayout.radiate(model)

        table = Dssr2DLayout._pair_table(model)
        try:
            points = DssrNaview().coordinates(table)
        except Exception as error:
            try:
                model.warnings.append("NAView fallback: %s" % str(error))
            except Exception:
                pass
            model._dssr2d_layout_variant = "radiate fallback"
            return Dssr2DLayout.radiate(model)

        scale = 1.65
        points = [(scale * point[0], scale * point[1]) for point in points]
        points, is_trna = Dssr2DLayout._standardize_trna_orientation(model, points)
        model._dssr2d_layout_variant = (
            "NAView standard tRNA cloverleaf" if is_trna else "NAView standard"
        )
        return Dssr2DLayout._center(points)

    @staticmethod
    def compute(model, algorithm):
        name = str(algorithm or "standard").strip().lower()
        if name in (
            "standard",
            "naview",
            "varna",
            "classic",
            "publication",
            "smart",
            "auto",
        ):
            return Dssr2DLayout._naview_layout(model)
        if name in ("legacy radiate", "legacy", "radiate", "radial"):
            model._dssr2d_layout_variant = "legacy radiate"
            return Dssr2DLayout.radiate(model)
        if name in ("circular", "circle"):
            model._dssr2d_layout_variant = "circular"
            return Dssr2DLayout.circular(model)
        if name in ("linear", "line"):
            model._dssr2d_layout_variant = "linear"
            return Dssr2DLayout.linear(model)
        return Dssr2DLayout._naview_layout(model)

    @staticmethod
    def radiate(model):
        topology = Dssr2DLayout._detect_trna_topology(model)
        if topology is not None:
            return Dssr2DLayout._trna_cloverleaf(model, topology)
        model._dssr2d_layout_variant = "general radiate"
        return Dssr2DLayout._radiate_general(model)

    smart = radiate

    @staticmethod
    def naview(model):
        return Dssr2DLayout._naview_layout(model)

    standard = naview


# ---------------------------------------------------------------------------
# 2D Graphics Items: Edges, PyMOL Color Palette, and Nucleotide Nodes
# ---------------------------------------------------------------------------


class Dssr2DEdgeItem(QtWidgets.QGraphicsPathItem):
    def __init__(
        self, node_a, node_b, kind="backbone", layer=0, lw="", linear_layout=False
    ):
        super().__init__()
        self.node_a, self.node_b = node_a, node_b
        self.kind, self.layer = str(kind), int(layer)
        self.lw, self.linear_layout = str(lw or ""), bool(linear_layout)
        self.setZValue(-5.0 if self.kind == "backbone" else -3.0)
        self._set_style()
        self.update_geometry()

    def _set_style(self):
        is_dark = getattr(getattr(self.node_a, "viewer", None), "is_dark", False)
        if self.kind == "backbone":
            # Slate gray backbone
            color = (
                QtGui.QColor(100, 116, 139) if is_dark else QtGui.QColor(148, 163, 184)
            )
            width = 1.0
            style = QtCore.Qt.SolidLine
        elif self.kind == "tertiary":
            # Tertiary pairs: luminous violet
            color = (
                QtGui.QColor(192, 132, 252) if is_dark else QtGui.QColor(124, 58, 237)
            )
            width = 1.35
            style = QtCore.Qt.DashLine
        elif self.layer > 0:
            # Pseudoknots: bright purple
            color = (
                QtGui.QColor(216, 180, 254) if is_dark else QtGui.QColor(147, 51, 234)
            )
            width = 1.45
            style = QtCore.Qt.DashLine
        else:
            # Canonical Watson-Crick/Wobble rungs: electric sky-blue in dark mode, royal cobalt in light mode
            color = QtGui.QColor(56, 189, 248) if is_dark else QtGui.QColor(29, 78, 216)
            width = 2.2
            style = QtCore.Qt.SolidLine

        self.setPen(
            QtGui.QPen(color, width, style, QtCore.Qt.RoundCap, QtCore.Qt.RoundJoin)
        )
        if self.lw:
            self.setToolTip("Base pair: %s" % self.lw)

    def _pen(self, color, width):
        style = QtCore.Qt.DashLine if self.kind == "tertiary" else QtCore.Qt.SolidLine
        return QtGui.QPen(color, width, style, QtCore.Qt.RoundCap, QtCore.Qt.RoundJoin)

    def _draw_paths(self, painter, pens):
        """Share painter setup across plain lines and gel passes."""
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        painter.setBrush(QtCore.Qt.NoBrush)
        for pen in pens:
            painter.setPen(pen)
            painter.drawPath(self.path())

    def paint(self, painter, option, widget=None):
        """Paint lines for backbone and base pairs."""
        try:
            gel = self.node_a.viewer.gel_style_enabled()
        except Exception:
            gel = False

        if not gel:
            return self._paint_flat(painter, option, widget)

        # Multi-pass gel mode line rendering
        if self.kind == "backbone":
            rgba, width = (112, 154, 186, 225), 1.65
        elif self.kind == "tertiary":
            rgba, width = (235, 111, 255, 220), 1.55
        elif self.layer > 0:
            palette = ((183, 123, 255, 240), (255, 112, 176, 240), (255, 190, 82, 240))
            rgba, width = palette[(self.layer - 1) % len(palette)], 2.05
        else:
            rgba, width = (78, 200, 255, 245), 2.15

        color = QtGui.QColor(*rgba)
        saved = False
        try:
            painter.save()
            saved = True
            glow = QtGui.QColor(color)
            glow.setAlpha(48)
            self._draw_paths(
                painter, (self._pen(glow, width + 4.2), self._pen(color, width))
            )
            painter.restore()
        except Exception:
            if saved:
                try:
                    painter.restore()
                except Exception:
                    pass
            return self._paint_flat(painter, option, widget)

    def _paint_flat(self, painter, option, widget=None):
        try:
            painter.save()
            self._draw_paths(painter, (self.pen(),))
            painter.restore()
        except Exception:
            try:
                painter.restore()
            except Exception:
                pass
            super().paint(painter, option, widget)

    def update_geometry(self):
        first = self.node_a.pos()
        second = self.node_b.pos()
        path = QtGui.QPainterPath()
        path.moveTo(first)

        if self.linear_layout and self.kind != "backbone":
            span = abs(int(self.node_a.nt_index) - int(self.node_b.nt_index))
            height = min(330.0, 25.0 + 6.0 * span)
            sign = -1.0 if int(self.layer) % 2 else 1.0
            if self.kind == "tertiary":
                sign *= -1.0
            path.quadTo(
                QtCore.QPointF(0.5 * (first.x() + second.x()), sign * height),
                second,
            )
        elif self.kind == "tertiary" or self.layer > 0:
            # Arc tertiary and pseudoknot pairs gracefully above the intervening structure
            delta_x = second.x() - first.x()
            delta_y = second.y() - first.y()
            distance = math.hypot(delta_x, delta_y)
            if distance <= 1.0e-9:
                path.lineTo(second)
            else:
                normal_x = -delta_y / distance
                normal_y = delta_x / distance
                curvature = min(80.0, max(25.0, 0.16 * distance))
                middle_x = 0.5 * (first.x() + second.x())
                middle_y = 0.5 * (first.y() + second.y())
                path.quadTo(
                    QtCore.QPointF(
                        middle_x + normal_x * curvature,
                        middle_y + normal_y * curvature,
                    ),
                    second,
                )
        else:
            # Straight lines for backbone and standard stem base pairs
            path.lineTo(second)

        self.setPath(path)


class Dssr2DNodeItem(QtWidgets.QGraphicsEllipseItem):
    RADIUS = 13.5  # 27 px diameter fits NAView 44 px spacing cleanly

    def __init__(self, viewer, nt, x, y):
        radius = self.RADIUS
        super().__init__(-radius, -radius, 2.0 * radius, 2.0 * radius)
        self.viewer, self.nt = viewer, nt
        self.nt_index = int(nt.get("index", 0))
        self.edge_items = []
        self._dragging = False
        self._drag_origin = None
        self._drag_starts = {}
        self._drag_before = None
        self._hover = False
        self._pressed = False
        self._visual_scale = 1.0
        self._scale_target = 1.0
        self._scale_velocity = 0.0
        self._drag_weights = {}
        self._last_drag_delta = QtCore.QPointF(0.0, 0.0)
        self._last_move_pos = None
        self._last_move_time = None
        self._drag_speed = 0.0
        self.setFlags(
            QtWidgets.QGraphicsItem.ItemIsSelectable
            | QtWidgets.QGraphicsItem.ItemIsFocusable
            | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges
        )
        self.setAcceptHoverEvents(True)
        self.setCursor(QtCore.Qt.ArrowCursor)
        self.setCacheMode(QtWidgets.QGraphicsItem.DeviceCoordinateCache)
        self.setPos(float(x), float(y))
        self.setZValue(5.0)
        self._apply_style(False)
        self._add_text()
        DssrUI.no_mouse(self.base_text_item)
        self.setToolTip(self._tooltip())

    def _base_fill(self):
        if self.viewer.gel_style_enabled():
            return QtGui.QColor(236, 245, 252)
        return DssrUtils.base_style(self.nt.get("base", ""), self.viewer.base_colors)[
            "fill"
        ]

    def _apply_style(self, selected):
        color = QtGui.QColor(225, 29, 72) if selected else QtGui.QColor(45, 45, 45)
        self.setPen(QtGui.QPen(color, 2.8 if selected else 1.2))
        self.setBrush(QtGui.QBrush(self._base_fill()))

    def _add_text(self):
        text = QtWidgets.QGraphicsSimpleTextItem(str(self.nt.get("base", "N")), self)
        font = QtGui.QFont("Sans Serif")
        font.setPointSize(12)  # High-legibility 12pt Bold
        font.setBold(True)
        text.setFont(font)
        rect = text.boundingRect()
        text.setPos(-rect.width() / 2.0, -rect.height() / 2.0)
        text.setBrush(
            QtGui.QBrush(
                DssrUtils.base_text_color(
                    self.nt.get("base", ""), self.viewer.base_colors
                )
            )
        )
        self.base_text_item = text

    def _tooltip(self):
        pieces = ["nt %d" % int(self.nt.get("number", self.nt_index + 1))]
        if self.nt.get("nt_id"):
            pieces.append(str(self.nt.get("nt_id")))
        else:
            if self.nt.get("chain"):
                pieces.append("chain %s" % self.nt.get("chain"))
            if self.nt.get("resi"):
                pieces.append("resi %s" % self.nt.get("resi"))
        return "\n".join(pieces)

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
            for edge in self.edge_items:
                edge.update_geometry()
        result = super().itemChange(change, value)
        try:
            if change == QtWidgets.QGraphicsItem.ItemSelectedHasChanged:
                selected = bool(value)
                self._apply_style(selected)
                self.setZValue(12.0 if selected else (16.0 if self._hover else 5.0))
                self._set_target_scale(
                    1.045 if selected else (1.075 if self._hover else 1.0),
                    kick=0.010 if selected else 0.0,
                )
                self.update()
            elif change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
                self.viewer._schedule_scene_rect()
        except Exception:
            pass
        return result

    def mousePressEvent(self, event):
        self._pressed = True
        self._set_target_scale(0.935, kick=-0.035)
        self._last_move_pos = event.scenePos()
        self._last_move_time = time.monotonic()
        self._drag_speed = 0.0
        button = event.button()
        if button != QtCore.Qt.LeftButton:
            self.viewer.select_nucleotide(self.nt_index)
            super().mousePressEvent(event)
            return

        modifiers = event.modifiers()
        if modifiers & QtCore.Qt.ControlModifier:
            self.setSelected(not self.isSelected())
            if not self.isSelected():
                self.viewer._sync_pymol_selection()
                event.accept()
                return
        elif modifiers & QtCore.Qt.ShiftModifier:
            self.setSelected(True)
        elif not self.isSelected():
            try:
                self.scene().clearSelection()
            except Exception:
                pass
            self.setSelected(True)
        selected = [node for node in self.viewer.nodes if node.isSelected()]
        if self not in selected:
            selected.append(self)
            self.setSelected(True)
        self._dragging = True
        self._drag_origin = event.scenePos()
        self._drag_before = self.viewer._capture_positions()
        self._drag_starts = {
            node.nt_index: QtCore.QPointF(node.pos()) for node in selected
        }
        try:
            self.setCursor(QtCore.Qt.ArrowCursor)
            self.viewer.view.setFocus()
        except Exception:
            pass
        self.viewer._sync_pymol_selection()
        event.accept()
        try:
            self.viewer._prepare_node_drag(self, modifiers)
        except Exception:
            pass

    def mouseMoveEvent(self, event):
        if not self._dragging or self._drag_origin is None:
            super().mouseMoveEvent(event)
            return
        delta = event.scenePos() - self._drag_origin
        self._last_drag_delta = QtCore.QPointF(delta)
        now = time.monotonic()
        if self._last_move_pos is not None and self._last_move_time is not None:
            dt = max(0.001, now - self._last_move_time)
            step = event.scenePos() - self._last_move_pos
            instant = math.hypot(step.x(), step.y()) / dt
            self._drag_speed = 0.72 * self._drag_speed + 0.28 * instant
        self._last_move_pos, self._last_move_time = event.scenePos(), now
        self._apply_drag_delta(delta, smooth=self.viewer.gel_style_enabled())
        event.accept()

    def _apply_drag_delta(self, delta, smooth=False):
        for index, start in self._drag_starts.items():
            if not 0 <= index < len(self.viewer.nodes):
                continue
            weight = float(self._drag_weights.get(index, 1.0))
            x = start.x() + delta.x() * weight
            y = start.y() + delta.y() * weight
            node = self.viewer.nodes[index]
            if smooth and weight < 0.999:
                follow = 0.44 + 0.34 * weight
                current = node.pos()
                x = current.x() + (x - current.x()) * follow
                y = current.y() + (y - current.y()) * follow
            node.setPos(x, y)

    def mouseReleaseEvent(self, event):
        was_dragging = self._dragging
        if was_dragging:
            self._apply_drag_delta(QtCore.QPointF(self._last_drag_delta))
        if self._dragging:
            self._dragging = False
            self.viewer._push_history(
                self._drag_before,
                self.viewer._capture_positions(),
                "move base%s" % ("s" if len(self._drag_starts) != 1 else ""),
            )
            self._drag_origin = self._drag_before = None
            self._drag_starts = {}
            try:
                self.setCursor(QtCore.Qt.ArrowCursor)
            except Exception:
                pass
            self.viewer._update_editor_status("manual move")
            event.accept()
        else:
            super().mouseReleaseEvent(event)
        self._pressed = False
        target = 1.075 if self._hover else (1.045 if self.isSelected() else 1.0)
        kick = (
            min(0.11, max(0.02, self._drag_speed * 0.00012)) if was_dragging else 0.025
        )
        self._set_target_scale(target, kick=kick)
        self._drag_weights = {}
        self._last_move_pos = self._last_move_time = None
        self.update()

    def contextMenuEvent(self, event):
        menu = QtWidgets.QMenu()
        header = menu.addAction(
            "%s  ·  nt %s"
            % (self.nt.get("base", "N"), self.nt.get("resi") or self.nt_index + 1)
        )
        header.setEnabled(False)
        menu.addSeparator()
        groups = {}
        for text, group in (
            ("Select base", lambda i: [i]),
            ("Select base pair", self.viewer._pair_group),
            ("Select loop / unpaired region", self.viewer._loop_group),
            ("Select stem", self.viewer._stem_indices),
            ("Select whole branch", self.viewer._branch_group),
        ):
            groups[menu.addAction(text)] = group
        menu.addSeparator()
        actions = {
            menu.addAction(text): callback
            for text, callback in (
                ("Center selection", self.viewer.fit_selected),
                (
                    "Reset selection to automatic layout",
                    self.viewer.reset_selected_bases,
                ),
                ("Undo", self.viewer.undo_layout),
            )
        }
        try:
            chosen = menu.exec_(event.screenPos())
        except Exception:
            try:
                chosen = menu.exec(event.screenPos())
            except Exception:
                chosen = None
        if chosen in groups:
            self.viewer.select_indices(groups[chosen](self.nt_index), replace=True)
        elif chosen in actions:
            actions[chosen]()
        event.accept()

    def boundingRect(self):
        radius = float(self.RADIUS)
        return QtCore.QRectF(
            -radius - 7.0,
            -radius - 7.0,
            2.0 * radius + 14.0,
            2.0 * radius + 14.0,
        )

    def shape(self):
        radius = float(self.RADIUS) + 3.0
        path = QtGui.QPainterPath()
        path.addEllipse(QtCore.QRectF(-radius, -radius, 2.0 * radius, 2.0 * radius))
        return path

    def paint(self, painter, option, widget=None):
        radius = float(self.RADIUS)
        gel = self.viewer.gel_style_enabled()
        saved = False
        try:
            painter.save()
            saved = True
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            selected = bool(self.isSelected())
            hovered = bool(getattr(self, "_hover", False))
            pressed = bool(getattr(self, "_pressed", False))

            # Drop shadow strictly in Gel mode; keep non-Gel completely flat and crisp
            if gel:
                shadow = QtGui.QColor(2, 8, 23, 105)
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(shadow))
                painter.drawEllipse(
                    QtCore.QRectF(
                        -radius + 2.2,
                        -radius + 3.4,
                        2.0 * radius,
                        2.0 * radius,
                    )
                )

            # Selection aura: PyMOL hot pink glow
            if selected or hovered:
                aura_color = (
                    QtGui.QColor(244, 63, 94, 125)
                    if selected
                    else QtGui.QColor(148, 163, 184, 70)
                )
                aura_radius = radius + (5.0 if selected else 3.0)
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(aura_color))
                painter.drawEllipse(
                    QtCore.QRectF(
                        -aura_radius,
                        -aura_radius,
                        2.0 * aura_radius,
                        2.0 * aura_radius,
                    )
                )

            show_circle = getattr(self.viewer, "show_circles", True)

            if gel:
                base = QtGui.QColor(self._base_fill())
                gradient = QtGui.QRadialGradient(
                    QtCore.QPointF(-radius * 0.38, -radius * 0.48),
                    radius * 1.58,
                )
                gradient.setColorAt(0.00, QtGui.QColor(255, 255, 255, 252))
                gradient.setColorAt(0.18, base.lighter(148))
                gradient.setColorAt(0.62, base.lighter(106))
                gradient.setColorAt(1.00, base.darker(132))
                fill = QtGui.QBrush(gradient)
                border = QtGui.QColor(188, 235, 255, 225)
            else:
                # Flat classical PyMOL base colors (soft pastel fill + crisp border)
                style = DssrUtils.base_style(
                    self.nt.get("base", ""), self.viewer.base_colors
                )
                fill = QtGui.QBrush(style["fill"])
                border = style["border"]

            if selected:
                # Signature PyMOL selection pink/crimson (#e11d48) with soft rose fill
                border = QtGui.QColor(225, 29, 72, 255)
                fill = QtGui.QBrush(QtGui.QColor(255, 228, 230))
            elif hovered:
                border = (
                    border.lighter(125)
                    if self.viewer.base_colors
                    else QtGui.QColor(105, 225, 255, 245)
                )
            if pressed:
                border = QtGui.QColor(15, 23, 42, 255)

            # Draw the circle if circles are enabled, or if the node is currently selected/hovered
            if show_circle or selected or hovered or gel:
                pen = QtGui.QPen(border)
                pen.setWidthF(2.4 if selected else (1.8 if hovered else 1.4))
                painter.setPen(pen)
                painter.setBrush(
                    fill if (show_circle or selected or gel) else QtCore.Qt.NoBrush
                )
                painter.drawEllipse(
                    QtCore.QRectF(-radius, -radius, 2.0 * radius, 2.0 * radius)
                )

            if gel:
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 110)))
                painter.drawEllipse(
                    QtCore.QRectF(
                        -radius * 0.58,
                        -radius * 0.68,
                        radius * 0.82,
                        radius * 0.38,
                    )
                )
            painter.restore()
        except Exception:
            if saved:
                try:
                    painter.restore()
                except Exception:
                    pass
            QtWidgets.QGraphicsEllipseItem.paint(self, painter, option, widget)

    def _set_target_scale(self, target, kick=0.0):
        self._scale_target = float(target)
        self._scale_velocity += float(kick)
        try:
            self.viewer._ensure_animation()
        except Exception:
            pass

    def _advance_visual(self):
        enabled = self.viewer.gel_style_enabled()
        target = self._scale_target if enabled else 1.0
        stiffness = 0.24 if enabled else 0.42
        damping = 0.68 if enabled else 0.55
        self._scale_velocity = (
            self._scale_velocity + (target - self._visual_scale) * stiffness
        ) * damping
        self._visual_scale += self._scale_velocity
        if (
            abs(target - self._visual_scale) < 0.0006
            and abs(self._scale_velocity) < 0.0006
        ):
            self._visual_scale = target
            self._scale_velocity = 0.0
        try:
            self.setScale(max(0.82, min(1.24, self._visual_scale)))
            self.update()
        except Exception:
            pass
        return not (
            abs(target - self._visual_scale) < 0.0007
            and abs(self._scale_velocity) < 0.0007
        )

    def hoverEnterEvent(self, event):
        self._hover = True
        self.setZValue(16.0)
        self._set_target_scale(1.075, kick=0.018)
        self.update()
        try:
            super().hoverEnterEvent(event)
        except Exception:
            pass

    def hoverLeaveEvent(self, event):
        self._hover = False
        self.setZValue(12.0 if self.isSelected() else 5.0)
        self._set_target_scale(1.045 if self.isSelected() else 1.0)
        self.update()
        try:
            super().hoverLeaveEvent(event)
        except Exception:
            pass


class Dssr2DGraphicsView(QtWidgets.QGraphicsView):
    """Empty-canvas rectangle selection, base editing, and keyboard-only pan."""

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.editor = parent
        self._rectangle_item = self._rectangle_origin = None
        self._rectangle_base = set()
        self._rectangle_mode = "replace"
        self._brushing = self._brush_erase = False
        self._brush_last = self._brush_path = self._brush_item = None
        self.setRenderHints(
            QtGui.QPainter.Antialiasing
            | QtGui.QPainter.TextAntialiasing
            | QtGui.QPainter.SmoothPixmapTransform
        )
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("white")))
        self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
        self.setInteractive(True)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorViewCenter)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.BoundingRectViewportUpdate)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setCursor(QtCore.Qt.ArrowCursor)

    def mousePressEvent(self, event):
        self.setFocus(QtCore.Qt.MouseFocusReason)
        if event.button() == QtCore.Qt.LeftButton:
            if self._tool() == "brush":
                self._begin_brush(event)
            elif self._node_item(self.itemAt(event.pos())) is None:
                self._begin_rectangle(event)
            else:
                super().mousePressEvent(event)
        elif event.button() == QtCore.Qt.MiddleButton:
            event.accept()
        else:
            # The normal scene route preserves a nucleotide's context menu.
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._rectangle_item is not None:
            self._update_rectangle(self.mapToScene(event.pos()))
            event.accept()
        elif self._brushing and self._brush_last is not None:
            point = self.mapToScene(event.pos())
            first, self._brush_last = self._brush_last, QtCore.QPointF(point)
            self._brush_path.lineTo(point)
            self._brush_item.setPath(self._brush_path)
            self._apply_brush_segment(first, point)
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton and self._rectangle_item is not None:
            self._update_rectangle(self.mapToScene(event.pos()))
            self._finish_rectangle()
            event.accept()
        elif event.button() == QtCore.Qt.LeftButton and self._brushing:
            self._finish_brush(event)
        elif event.button() == QtCore.Qt.MiddleButton:
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        key, modifiers = event.key(), event.modifiers()
        control = bool(modifiers & QtCore.Qt.ControlModifier)
        shift = bool(modifiers & QtCore.Qt.ShiftModifier)
        modes = dict(
            zip(
                (
                    QtCore.Qt.Key_1,
                    QtCore.Qt.Key_2,
                    QtCore.Qt.Key_3,
                    QtCore.Qt.Key_4,
                    QtCore.Qt.Key_5,
                    QtCore.Qt.Key_6,
                ),
                ("base", "selection", "pair", "loop", "stem", "branch"),
            )
        )
        arrows = {
            QtCore.Qt.Key_Left: (-1, 0),
            QtCore.Qt.Key_Right: (1, 0),
            QtCore.Qt.Key_Up: (0, -1),
            QtCore.Qt.Key_Down: (0, 1),
        }
        movement = dict(arrows)
        movement.update(
            {
                QtCore.Qt.Key_A: (-1, 0),
                QtCore.Qt.Key_D: (1, 0),
                QtCore.Qt.Key_W: (0, -1),
                QtCore.Qt.Key_S: (0, 1),
            }
        )
        if control and key == QtCore.Qt.Key_Z:
            (self.editor.redo_layout if shift else self.editor.undo_layout)()
        elif control and key == QtCore.Qt.Key_Y:
            self.editor.redo_layout()
        elif control and key == QtCore.Qt.Key_A:
            self.editor.select_all_bases()
        elif control and key in arrows:
            dx, dy = arrows[key]
            step = 10.0 if shift else 2.0
            self.editor.nudge_selected(dx * step, dy * step)
        elif key in movement and not control and not modifiers & QtCore.Qt.AltModifier:
            dx, dy = movement[key]
            step = 120 if shift else 40
            self._pan_view(dx * step, dy * step)
        elif key in (QtCore.Qt.Key_B, QtCore.Qt.Key_P):
            self.cancel_selection_gesture()
            self.editor.set_interaction_tool(
                "brush" if key == QtCore.Qt.Key_B else "edit"
            )
        elif key in modes:
            self.editor.set_drag_mode(modes[key])
        elif key == QtCore.Qt.Key_F:
            self.editor.fit_scene()
        elif key == QtCore.Qt.Key_C:
            self.editor.fit_selected()
        elif key == QtCore.Qt.Key_Escape:
            if self._rectangle_item is not None:
                self._finish_rectangle(cancel=True)
            else:
                self.cancel_selection_gesture()
                self.editor.clear_base_selection()
        elif key != QtCore.Qt.Key_Space:
            super().keyPressEvent(event)
            return
        event.accept()

    def _pan_view(self, dx, dy):
        """Move the drawing by screen pixels, including immediately after Fit."""
        scale = max(0.08, abs(self.transform().m11()))
        offset = QtCore.QPointF(dx / scale, dy / scale)
        anchor = self.mapToScene(QtCore.QPoint(0, 0))
        visible = self.mapToScene(self.viewport().rect()).boundingRect()
        margin = 100.0 / scale
        needed = visible.translated(-offset).adjusted(-margin, -margin, margin, margin)
        self.setSceneRect(self.sceneRect().united(needed))
        mapped = self.mapFromScene(anchor)
        horizontal, vertical = self.horizontalScrollBar(), self.verticalScrollBar()
        horizontal.setValue(horizontal.value() + mapped.x() - int(dx))
        vertical.setValue(vertical.value() + mapped.y() - int(dy))

    def _begin_rectangle(self, event):
        self.cancel_selection_gesture()
        self._rectangle_origin = self.mapToScene(event.pos())
        self._rectangle_base = {
            node.nt_index for node in self.editor.nodes if node.isSelected()
        }
        modifiers = event.modifiers()
        self._rectangle_mode = (
            "subtract"
            if modifiers & QtCore.Qt.AltModifier
            else "add"
            if modifiers & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
            else "replace"
        )
        item = QtWidgets.QGraphicsRectItem()
        pen = QtGui.QPen(QtGui.QColor(45, 126, 225, 230), 1.25)
        pen.setCosmetic(True)
        item.setPen(pen)
        item.setBrush(QtGui.QBrush(QtGui.QColor(66, 153, 225, 45)))
        item.setZValue(1000.0)
        item.setAcceptedMouseButtons(QtCore.Qt.NoButton)
        self.scene().addItem(item)
        self._rectangle_item = item
        self._update_rectangle(self._rectangle_origin)
        event.accept()

    def _set_rectangle_selection(self, wanted):
        changed = False
        was_rebuilding = self.editor._rebuilding
        self.editor._rebuilding = True
        try:
            for node in self.editor.nodes:
                selected = node.nt_index in wanted
                if node.isSelected() != selected:
                    node.setSelected(selected)
                    changed = True
        finally:
            self.editor._rebuilding = was_rebuilding
        # Throttle instead of restarting the debounce on every mouse move:
        # continuous rectangle motion still produces live 3D updates.
        if changed and not self.editor._sync_timer.isActive():
            self.editor._schedule_live_sync()

    def _update_rectangle(self, point):
        rect = QtCore.QRectF(self._rectangle_origin, point).normalized()
        self._rectangle_item.setRect(rect)
        inside = {
            node.nt_index
            for node in self.editor.nodes
            if rect.contains(node.scenePos())
        }
        if self._rectangle_mode == "add":
            wanted = self._rectangle_base | inside
        elif self._rectangle_mode == "subtract":
            wanted = self._rectangle_base - inside
        else:
            wanted = inside
        self._set_rectangle_selection(wanted)
        self.editor._update_editor_status("rectangle: %d bases" % len(wanted))

    def _finish_rectangle(self, cancel=False):
        if cancel:
            self._set_rectangle_selection(self._rectangle_base)
        self.cancel_selection_gesture()
        self.editor._flush_live_sync(final=True)
        self.editor._update_editor_status(
            "rectangle canceled" if cancel else "rectangle selection mapped to 3D"
        )

    def cancel_selection_gesture(self):
        """Drop temporary overlays/state without scheduling hidden-view work."""
        self.editor._sync_timer.stop()
        self.editor._sync_pending = False
        for item in (self._rectangle_item, self._brush_item):
            if item is not None:
                try:
                    if item.scene() is self.scene():
                        self.scene().removeItem(item)
                except RuntimeError:
                    pass  # The scene may already have disposed its old items.
        self._rectangle_item = self._rectangle_origin = None
        self._rectangle_base = set()
        self._rectangle_mode = "replace"
        self._brushing = self._brush_erase = False
        self._brush_last = self._brush_path = self._brush_item = None
        self.setCursor(
            QtCore.Qt.CrossCursor if self._tool() == "brush" else QtCore.Qt.ArrowCursor
        )

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        if delta:
            current = abs(float(self.transform().m11()))
            target = max(0.08, min(14.0, current * math.pow(1.00105, delta)))
            factor = target / max(1.0e-9, current)
            self.scale(factor, factor)
            event.accept()

    def mouseDoubleClickEvent(self, event):
        node = self._node_item(self.itemAt(event.pos()))
        if node is not None:
            self.editor.select_indices(
                self.editor._stem_indices(node.nt_index), replace=True
            )
            self.editor.fit_selected()
        else:
            self.editor.fit_scene()
        event.accept()

    @staticmethod
    def _node_item(item):
        for _ in range(6):
            if item is None or isinstance(item, Dssr2DNodeItem):
                return item
            item = item.parentItem()
        return None

    def _tool(self):
        return self.editor.interaction_tool()

    def _scene_radius(self):
        return self.editor.brush_radius_spin.value() / max(
            0.08, abs(self.transform().m11())
        )

    @staticmethod
    def _distance_to_segment(point, first, second):
        dx, dy = second.x() - first.x(), second.y() - first.y()
        length_sq = dx * dx + dy * dy
        if length_sq <= 1.0e-12:
            return math.hypot(point.x() - first.x(), point.y() - first.y())
        ratio = (
            (point.x() - first.x()) * dx + (point.y() - first.y()) * dy
        ) / length_sq
        ratio = max(0.0, min(1.0, ratio))
        return math.hypot(
            point.x() - (first.x() + ratio * dx), point.y() - (first.y() + ratio * dy)
        )

    def _begin_brush(self, event):
        self.cancel_selection_gesture()
        self._brushing = True
        self._brush_erase = bool(event.modifiers() & QtCore.Qt.AltModifier)
        additive = bool(
            event.modifiers() & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
        )
        point = self.mapToScene(event.pos())
        self._brush_last = QtCore.QPointF(point)
        if not additive and not self._brush_erase:
            self.editor._rebuilding = True
            try:
                self.scene().clearSelection()
            finally:
                self.editor._rebuilding = False
        self._brush_path = QtGui.QPainterPath(point)
        item = QtWidgets.QGraphicsPathItem(self._brush_path)
        color = (
            QtGui.QColor(255, 100, 176, 70)
            if self._brush_erase
            else QtGui.QColor(58, 214, 255, 68)
        )
        pen = QtGui.QPen(color)
        pen.setWidthF(2.0 * self._scene_radius())
        pen.setCapStyle(QtCore.Qt.RoundCap)
        pen.setJoinStyle(QtCore.Qt.RoundJoin)
        item.setPen(pen)
        item.setBrush(QtGui.QBrush(QtCore.Qt.NoBrush))
        item.setZValue(3.0)
        item.setAcceptedMouseButtons(QtCore.Qt.NoButton)
        self.scene().addItem(item)
        self._brush_item = item
        self._apply_brush_segment(point, point)
        self.setCursor(QtCore.Qt.CrossCursor)
        event.accept()

    def _apply_brush_segment(self, first, second):
        radius, changed = self._scene_radius(), False
        self.editor._rebuilding = True
        try:
            for node in self.editor.nodes:
                if (
                    self._distance_to_segment(node.scenePos(), first, second)
                    <= radius + node.RADIUS * 0.65
                ):
                    wanted = not self._brush_erase
                    if node.isSelected() != wanted:
                        node.setSelected(wanted)
                        changed = True
        finally:
            self.editor._rebuilding = False
        if changed:
            self.editor._after_brush_selection(final=False)

    def _finish_brush(self, event=None):
        if not self._brushing:
            return
        self._brushing, self._brush_last = False, None
        if self._brush_item is not None:
            self.scene().removeItem(self._brush_item)
        self._brush_item = self._brush_path = None
        self.editor._after_brush_selection(final=True)
        self.editor._update_view_cursor()
        if event is not None:
            event.accept()


class Dssr2DSequenceView(QtWidgets.QTextEdit):
    """One text document for the sequence, real residue ruler, and linked selection."""

    def __init__(self, editor):
        super().__init__(editor)
        self.editor = editor
        self._selected = frozenset()
        self._anchor = self._drag_anchor = None
        self._gesture_base = set()
        self._gesture_toggle = False
        self._letter_colors = None
        self.setReadOnly(True)
        self.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setFocusPolicy(QtCore.Qt.NoFocus)
        self.setMouseTracking(True)
        self.viewport().setCursor(QtCore.Qt.ArrowCursor)
        font = QtGui.QFont()
        font.setStyleHint(QtGui.QFont.Monospace)
        font.setPointSize(10)
        font.setFamilies(["Menlo", "Monaco", "Courier New", "DejaVu Sans Mono"])
        self.setFont(font)
        self.document().setDocumentMargin(8)
        self.setFixedHeight(
            self.fontMetrics().height() * 2
            + 24
            + self.style().pixelMetric(QtWidgets.QStyle.PM_ScrollBarExtent)
        )
        self.setStyleSheet(
            "QTextEdit { background: white; color: #64748b; border: 1px solid #dce5ea; border-radius: 5px; }"
        )
        self.setToolTip(
            "Sequence: click or drag a range · Shift: extend · Ctrl: toggle\nNumbers are PyMOL residue IDs; hover for the original DSSR identifier."
        )
        self._build_document()
        self.set_letter_colors(editor.base_colors)

    @staticmethod
    def _text_length(text):
        return len(text.encode("utf-16-le")) // 2

    def _build_document(self):
        parts, spans, rulers = [], [], []
        character_offset = document_offset = 0
        nts, breaks = self.editor.model.nts, self.editor.model.chain_breaks
        for index, nt in enumerate(nts):
            chain_start = index == 0 or index - 1 in breaks
            prefix = (
                (("  |  " if index else "") + (str(nt.get("chain", "")) or "–") + ": ")
                if chain_start
                else " "
            )
            parts.append(prefix)
            character_offset += len(prefix)
            document_offset += self._text_length(prefix)
            base = str(nt.get("base", ""))
            start = document_offset
            parts.append(base)
            document_offset += self._text_length(base)
            spans.append((start, document_offset))
            if chain_start or (index + 1) % 10 == 0 or index == len(nts) - 1:
                residue = str(nt.get("resi", "")).strip()
                if residue:
                    rulers.append((character_offset, residue))
            character_offset += len(base)
        sequence = "".join(parts)
        ruler = [" "] * len(sequence)
        occupied = -1
        for start, label in rulers:
            if start <= occupied:
                continue
            end = start + len(label)
            if end > len(ruler):
                ruler.extend(" " * (end - len(ruler)))
            ruler[start:end] = label
            occupied = end
        ruler = "".join(ruler)
        self.setPlainText(ruler + "\n" + sequence)
        offset = self._text_length(ruler) + 1
        self._spans = [(start + offset, end + offset) for start, end in spans]
        self._starts = [start for start, _end in self._spans]

    def set_letter_colors(self, enabled):
        enabled = bool(enabled)
        if enabled == self._letter_colors:
            return
        self._letter_colors = enabled
        cursor = QtGui.QTextCursor(self.document())
        cursor.beginEditBlock()
        for nt, (start, end) in zip(self.editor.model.nts, self._spans):
            cursor.setPosition(start)
            cursor.setPosition(end, QtGui.QTextCursor.KeepAnchor)
            style = QtGui.QTextCharFormat()
            style.setForeground(DssrUtils.base_text_color(nt.get("base", ""), enabled))
            style.setFontWeight(QtGui.QFont.Bold)
            cursor.setCharFormat(style)
        cursor.endEditBlock()

    def set_selected_indices(self, indices):
        wanted = frozenset(index for index in indices if 0 <= index < len(self._spans))
        if wanted == self._selected:
            return
        self._selected = wanted
        groups = []
        for index in sorted(wanted):
            if groups and index == groups[-1][1] + 1:
                groups[-1][1] = index
            else:
                groups.append([index, index])
        selections = []
        for first, last in groups:
            selection = QtWidgets.QTextEdit.ExtraSelection()
            selection.cursor = QtGui.QTextCursor(self.document())
            selection.cursor.setPosition(self._spans[first][0])
            selection.cursor.setPosition(
                self._spans[last][1], QtGui.QTextCursor.KeepAnchor
            )
            is_dark = getattr(self.editor, "is_dark", False)
            selection.format.setBackground(
                QtGui.QColor("#881337") if is_dark else QtGui.QColor("#fce7f3")
            )
            selections.append(selection)
        self.setExtraSelections(selections)

    def _index_at(self, position):
        offset = self.cursorForPosition(position).position()
        index = bisect.bisect_right(self._starts, offset) - 1
        return index if index >= 0 and offset <= self._spans[index][1] else None

    def _apply_range(self, index):
        indices = set(
            range(min(self._drag_anchor, index), max(self._drag_anchor, index) + 1)
        )
        wanted = (
            self._gesture_base.symmetric_difference(indices)
            if self._gesture_toggle
            else indices
        )
        if wanted != self._selected:
            self.editor.select_indices(wanted, replace=True)

    def mousePressEvent(self, event):
        index = self._index_at(event.pos())
        if event.button() != QtCore.Qt.LeftButton or index is None:
            super().mousePressEvent(event)
            return
        self.editor.view.setFocus(QtCore.Qt.MouseFocusReason)
        shift = bool(event.modifiers() & QtCore.Qt.ShiftModifier)
        self._drag_anchor = (
            self._anchor if shift and self._anchor is not None else index
        )
        if not shift or self._anchor is None:
            self._anchor = index
        self._gesture_base = set(self._selected)
        self._gesture_toggle = bool(event.modifiers() & QtCore.Qt.ControlModifier)
        self._apply_range(index)
        event.accept()

    def mouseMoveEvent(self, event):
        index = self._index_at(event.pos())
        if index is not None:
            nt = self.editor.model.nts[index]
            self.viewport().setToolTip(
                "%s · %s\n%s"
                % (
                    nt.get("base", ""),
                    nt.get("name", ""),
                    nt.get("nt_id")
                    or "Sequence position %d; no DSSR residue identifier" % (index + 1),
                )
            )
        if self._drag_anchor is not None:
            if index is not None:
                self._apply_range(index)
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton and self._drag_anchor is not None:
            index = self._index_at(event.pos())
            if index is not None:
                self._apply_range(index)
            self._drag_anchor = None
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def wheelEvent(self, event):
        delta = event.angleDelta().x() or event.angleDelta().y()
        bar = self.horizontalScrollBar()
        bar.setValue(bar.value() - delta)
        event.accept()


class Dssr2DEditor(QtWidgets.QWidget):
    HISTORY_LIMIT = 100

    HIGHLIGHT_OBJECT = "_dssr_2d_brush_highlight"

    def __init__(
        self,
        model,
        pymol_selection="all",
        algorithm="standard",
        number_every=10,
        show_tertiary=False,
        parent=None,
    ):
        super().__init__(parent)
        self.model = model
        self.pymol_selection = str(pymol_selection or "all")
        requested = str(algorithm or "standard").strip().lower()
        self.algorithm = requested if requested in LAYOUT_CHOICES else "standard"
        self.number_every = max(0, int(number_every))
        self.show_tertiary = bool(show_tertiary)
        self.base_colors = True
        self.show_circles = True  # Default: circles visible
        self.is_dark = False
        self.nodes, self.edges = [], []
        self._rebuilding = False
        self._auto_positions = []
        self._undo, self._redo = [], []
        self._scene_rect_pending = False
        self._pair_table = Dssr2DLayout._planar_pair_table(model)
        self._stems = Dssr2DLayout._stem_tree(self._pair_table, model.chain_breaks)
        self._adjacency = None
        self._sync_pending = False
        self._sync_from_pymol = False
        self._last_pymol_signature = None
        self._last_highlight_signature = None
        self._last_sel_count = -1
        self._last_active_names = ()
        self._closed = False
        self._view_active = True
        self._shown_once = False

        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self._animation_tick)
        self._sync_timer = QtCore.QTimer(self)
        self._sync_timer.setSingleShot(True)
        self._sync_timer.setInterval(55)
        self._sync_timer.timeout.connect(self._flush_live_sync)
        self._reverse_timer = QtCore.QTimer(self)
        self._reverse_timer.setInterval(450)
        self._reverse_timer.timeout.connect(self._pull_pymol_selection)

        # Build widgets once and establish initial theme
        self._build_widgets()
        self.set_theme(False)

        self.scene.selectionChanged.connect(self._selection_changed)
        self.layout_combo.currentTextChanged.connect(self.redraw)
        self.number_spin.valueChanged.connect(self.redraw)
        self.tertiary_cb.toggled.connect(self.redraw)
        self.base_colors_cb.toggled.connect(self.redraw)
        self.redraw()
        self._update_history_buttons()
        self._update_view_cursor()
        self._update_editor_status("ready")
        self._reverse_timer.start()
        QtCore.QTimer.singleShot(0, self.fit_scene)

    def set_theme(self, is_dark):
        self.is_dark = bool(is_dark)
        self.setStyleSheet(DARK_THEME if self.is_dark else LIGHT_THEME)
        bg = QtGui.QColor("#0f172a" if self.is_dark else "#ffffff")
        self.view.setBackgroundBrush(QtGui.QBrush(bg))
        seq_bg = "#0f172a" if self.is_dark else "#ffffff"
        seq_text = "#94a3b8" if self.is_dark else "#475569"
        seq_border = "#334155" if self.is_dark else "#cbd5e1"
        self.sequence_view.setStyleSheet(
            "QTextEdit { background: %s; color: %s; border: 1px solid %s; border-radius: 5px; }"
            % (seq_bg, seq_text, seq_border)
        )
        for edge in self.edges:
            edge._set_style()
        self._refresh_scene_style()

    @staticmethod
    def _segment_intersects_box(p1, p2, rx1, ry1, rx2, ry2):
        """Liang-Barsky line clipping: returns True if segment p1-p2 intersects box."""
        if rx1 > rx2:
            rx1, rx2 = rx2, rx1
        if ry1 > ry2:
            ry1, ry2 = ry2, ry1

        min_x, max_x = (p1[0], p2[0]) if p1[0] <= p2[0] else (p2[0], p1[0])
        min_y, max_y = (p1[1], p2[1]) if p1[1] <= p2[1] else (p2[1], p1[1])
        if max_x < rx1 or min_x > rx2 or max_y < ry1 or min_y > ry2:
            return False

        if rx1 <= p1[0] <= rx2 and ry1 <= p1[1] <= ry2:
            return True
        if rx1 <= p2[0] <= rx2 and ry1 <= p2[1] <= ry2:
            return True

        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        p = [-dx, dx, -dy, dy]
        q = [p1[0] - rx1, rx2 - p1[0], p1[1] - ry1, ry2 - p1[1]]

        u1, u2 = 0.0, 1.0
        for pi, qi in zip(p, q):
            if pi == 0:
                if qi < 0:
                    return False
            else:
                t = qi / pi
                if pi < 0:
                    if t > u2:
                        return False
                    if t > u1:
                        u1 = t
                else:
                    if t < u1:
                        return False
                    if t < u2:
                        u2 = t
        return u1 <= u2

    def _add_edge(self, i, j, kind, layer=0, lw="", linear=False):
        if i < 0 or j < 0 or i >= len(self.nodes) or j >= len(self.nodes):
            return None
        edge = Dssr2DEdgeItem(
            self.nodes[i],
            self.nodes[j],
            kind=kind,
            layer=layer,
            lw=lw,
            linear_layout=linear,
        )
        self.scene.addItem(edge)
        self.nodes[i].edge_items.append(edge)
        self.nodes[j].edge_items.append(edge)
        self.edges.append(edge)
        edge.setAcceptedMouseButtons(QtCore.Qt.NoButton)
        edge.setAcceptHoverEvents(False)
        return edge

    def redraw(self, *_args):
        if self._closed or self._rebuilding:
            return
        self.view.cancel_selection_gesture()
        self._rebuilding = True
        try:
            self.algorithm = (
                self.layout_combo.currentText().strip().lower() or "standard"
            )
            self.number_every = self.number_spin.value()
            self.show_tertiary = self.tertiary_cb.isChecked()
            self.base_colors = self.base_colors_cb.isChecked()
            positions = self._capture_positions()
            selected = [node.nt_index for node in self.nodes if node.isSelected()]
            preserve = bool(positions) and self.sender() not in (
                self.layout_combo,
                self.redraw_btn,
            )
            if preserve:
                transform = QtGui.QTransform(self.view.transform())
                center = self.view.mapToScene(self.view.viewport().rect().center())
            else:
                positions = Dssr2DLayout.compute(self.model, self.algorithm)
                self._auto_positions = [(float(x), float(y)) for x, y in positions]
            self.scene.blockSignals(True)
            self.scene.clear()
            self.nodes, self.edges = [], []
            for nt, (x, y) in zip(self.model.nts, positions):
                node = Dssr2DNodeItem(self, nt, x, y)
                self.scene.addItem(node)
                self.nodes.append(node)
            linear = self.algorithm == "linear"
            for index in range(len(self.nodes) - 1):
                if index not in self.model.chain_breaks:
                    self._add_edge(index, index + 1, "backbone", linear=linear)
            groups = [("secondary", self.model.secondary_pairs)]
            if self.show_tertiary:
                groups.append(("tertiary", self.model.tertiary_pairs))
            for kind, pairs in groups:
                for pair in pairs:
                    layer = -1 if kind == "tertiary" else int(pair.get("layer", 0))
                    self._add_edge(
                        int(pair["i"]),
                        int(pair["j"]),
                        kind,
                        layer=layer,
                        lw=pair.get("lw", ""),
                        linear=linear,
                    )
            self._chain_rects = self._add_chain_labels()
            self._add_number_labels()
            self._update_scene_rect()
            for index in selected:
                if index < len(self.nodes):
                    self.nodes[index].setSelected(True)
            self.scene.blockSignals(False)
            if preserve:
                self.view.setTransform(transform)
                self.view.centerOn(center)
            else:
                self.fit_scene()
            self._update_editor_status("redrawn" if preserve else "automatic layout")
        finally:
            self.scene.blockSignals(False)
            self._rebuilding = False
        for node in self.nodes:
            node.setScale(1.0)
        self._refresh_scene_style()
        self._ensure_animation()

    def _add_number_labels(self):
        """Add sparse residue numbers strictly outside helices and clear of all linkages."""
        total = len(self.nodes)
        if total <= 0:
            return
        period = int(self.number_every)
        indices = set()
        if period > 0:
            indices.update(index for index in range(total) if (index + 1) % period == 0)
        for break_after in self.model.chain_breaks:
            if 0 <= break_after < total:
                indices.add(break_after)
            if 0 <= break_after + 1 < total:
                indices.add(break_after + 1)

        center_x = sum(node.pos().x() for node in self.nodes) / float(total)
        center_y = sum(node.pos().y() for node in self.nodes) / float(total)
        node_radius = float(getattr(Dssr2DNodeItem, "RADIUS", 13.5)) + 2.5
        node_rects = [
            QtCore.QRectF(
                node.pos().x() - node_radius,
                node.pos().y() - node_radius,
                2.0 * node_radius,
                2.0 * node_radius,
            )
            for node in self.nodes
        ]

        # Collect all physical edge lines (backbone and rungs) to protect them from overlap
        edge_segments = []
        for i in range(total - 1):
            if i not in self.model.chain_breaks:
                p1 = self.nodes[i].pos()
                p2 = self.nodes[i + 1].pos()
                edge_segments.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

        for pair in self.model.secondary_pairs:
            i = int(pair.get("i", -1))
            j = int(pair.get("j", -1))
            if 0 <= i < total and 0 <= j < total and i != j:
                p1 = self.nodes[i].pos()
                p2 = self.nodes[j].pos()
                edge_segments.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

        if self.show_tertiary:
            for pair in self.model.tertiary_pairs:
                i = int(pair.get("i", -1))
                j = int(pair.get("j", -1))
                if 0 <= i < total and 0 <= j < total and i != j:
                    p1 = self.nodes[i].pos()
                    p2 = self.nodes[j].pos()
                    edge_segments.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

        used_label_rects = list(getattr(self, "_chain_rects", []))

        def _intersection_area(first, second):
            try:
                overlap = first.intersected(second)
                if overlap.isEmpty():
                    return 0.0
                return max(0.0, overlap.width()) * max(0.0, overlap.height())
            except Exception:
                return 0.0

        is_dark = getattr(self, "is_dark", False)
        num_color = QtGui.QColor(248, 250, 252) if is_dark else QtGui.QColor(15, 23, 42)

        for index in sorted(indices):
            node = self.nodes[index]
            nt = self.model.nts[index]
            value = str(nt.get("resi") or nt.get("number", index + 1))
            label = QtWidgets.QGraphicsSimpleTextItem(value, node)
            font = QtGui.QFont("Sans Serif")
            font.setPointSize(11)
            font.setBold(True)
            label.setFont(font)
            label.setBrush(QtGui.QBrush(num_color))

            pos = node.pos()
            partner_idx = self._pair_partner(index)

            # 1. Determine strict OUTWARD vector
            if 0 <= partner_idx < total:
                # Paired base: points strictly away from base-pair partner
                p_pos = self.nodes[partner_idx].pos()
                vx = pos.x() - p_pos.x()
                vy = pos.y() - p_pos.y()
                vlen = math.hypot(vx, vy)
                outward_x, outward_y = (
                    (vx / vlen, vy / vlen) if vlen > 1e-6 else (1.0, 0.0)
                )
            else:
                # Unpaired base (loops/linkers): use curve curvature (2N - P - F)
                previous = (
                    self.nodes[index - 1].pos()
                    if index > 0 and (index - 1) not in self.model.chain_breaks
                    else None
                )
                following = (
                    self.nodes[index + 1].pos()
                    if index + 1 < total and index not in self.model.chain_breaks
                    else None
                )

                if previous is not None and following is not None:
                    kx = 2.0 * pos.x() - previous.x() - following.x()
                    ky = 2.0 * pos.y() - previous.y() - following.y()
                    klen = math.hypot(kx, ky)
                    if klen > 0.5:
                        outward_x, outward_y = kx / klen, ky / klen
                    else:
                        tx = following.x() - previous.x()
                        ty = following.y() - previous.y()
                        nx, ny = -ty, tx
                        nlen = math.hypot(nx, ny)
                        if nlen > 1e-6:
                            nx /= nlen
                            ny /= nlen
                            if (
                                nx * (pos.x() - center_x) + ny * (pos.y() - center_y)
                                < 0.0
                            ):
                                nx, ny = -nx, -ny
                            outward_x, outward_y = nx, ny
                        else:
                            outward_x, outward_y = 1.0, 0.0
                elif previous is not None:
                    tx = pos.x() - previous.x()
                    ty = pos.y() - previous.y()
                    outward_x, outward_y = (
                        tx / max(1e-6, math.hypot(tx, ty)),
                        ty / max(1e-6, math.hypot(tx, ty)),
                    )
                elif following is not None:
                    tx = following.x() - pos.x()
                    ty = following.y() - pos.y()
                    outward_x, outward_y = (
                        -tx / max(1e-6, math.hypot(tx, ty)),
                        -ty / max(1e-6, math.hypot(tx, ty)),
                    )
                else:
                    cx = pos.x() - center_x
                    cy = pos.y() - center_y
                    clen = math.hypot(cx, cy)
                    outward_x, outward_y = (
                        (cx / clen, cy / clen) if clen > 1e-6 else (1.0, 0.0)
                    )

            # 2. Candidate directions: only outward hemisphere (within +/- 60 deg of outward)
            candidate_angles = [0.0, 0.35, -0.35, 0.70, -0.70, 1.05, -1.05]
            candidate_distances = (22.0, 28.0, 35.0, 42.0)
            rect = label.boundingRect()
            best = None

            for dir_rank, alpha in enumerate(candidate_angles):
                cosa = math.cos(alpha)
                sina = math.sin(alpha)
                dx = outward_x * cosa - outward_y * sina
                dy = outward_x * sina + outward_y * cosa

                for dist_rank, dist in enumerate(candidate_distances):
                    lx = dx * dist - 0.5 * rect.width()
                    ly = dy * dist - 0.5 * rect.height()
                    srect = QtCore.QRectF(
                        pos.x() + lx, pos.y() + ly, rect.width(), rect.height()
                    )

                    score = 0.05 * dist_rank + 0.02 * dir_rank

                    # Avoid other node bubbles
                    for n_idx, n_rect in enumerate(node_rects):
                        if n_idx == index:
                            continue
                        area = _intersection_area(srect, n_rect)
                        if area > 0:
                            score += 25.0 * area

                    # Avoid previously placed text
                    for u_rect in used_label_rects:
                        area = _intersection_area(srect, u_rect)
                        if area > 0:
                            score += 50.0 * area

                    # Avoid backbone and base-pair linkages (3 px buffer)
                    for p1, p2 in edge_segments:
                        if self._segment_intersects_box(
                            p1,
                            p2,
                            srect.left() - 3.0,
                            srect.top() - 3.0,
                            srect.right() + 3.0,
                            srect.bottom() + 3.0,
                        ):
                            score += 60.0

                    if best is None or score < best[0]:
                        best = (score, lx, ly, srect)

            if best is None:
                best = (0.0, outward_x * 24.0, outward_y * 24.0, QtCore.QRectF())

            label.setPos(best[1], best[2])
            used_label_rects.append(best[3])
            label.setZValue(8.0)
            DssrUI.no_mouse(label)

    def _add_chain_labels(self):
        total = len(self.nodes)
        if total <= 0:
            return []
        segments = []
        start = 0
        for break_after in sorted(self.model.chain_breaks):
            if break_after >= start:
                segments.append((start, min(total - 1, break_after)))
                start = break_after + 1
        if start < total:
            segments.append((start, total - 1))
        if not segments:
            segments = [(0, total - 1)]

        is_dark = getattr(self, "is_dark", False)
        term_color = QtGui.QColor(56, 189, 248) if is_dark else QtGui.QColor(15, 23, 42)
        chain_rects = []

        for segment_number, (first, last) in enumerate(segments, 1):
            for index, text_value, offset in (
                (first, "5′", (-42.0, -26.0)),
                (last, "3′", (24.0, -26.0)),
            ):
                label = QtWidgets.QGraphicsSimpleTextItem(text_value, self.nodes[index])
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(13)
                font.setBold(True)
                label.setFont(font)
                label.setBrush(QtGui.QBrush(term_color))
                label.setPos(offset[0], offset[1])
                label.setZValue(8.0)
                DssrUI.no_mouse(label)
                br = label.boundingRect()
                pos = self.nodes[index].pos()
                chain_rects.append(
                    QtCore.QRectF(
                        pos.x() + offset[0],
                        pos.y() + offset[1],
                        br.width(),
                        br.height(),
                    )
                )

            if len(segments) > 1:
                chain = str(self.model.nts[first].get("chain", "")).strip()
                text_value = "chain %s" % (chain or segment_number)
                label = QtWidgets.QGraphicsSimpleTextItem(text_value, self.nodes[first])
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(11)
                font.setBold(True)
                label.setFont(font)
                label.setBrush(QtGui.QBrush(term_color))
                label.setPos(-48.0, -52.0)
                label.setZValue(8.0)
                DssrUI.no_mouse(label)
                br = label.boundingRect()
                pos = self.nodes[first].pos()
                chain_rects.append(
                    QtCore.QRectF(
                        pos.x() - 48.0,
                        pos.y() - 52.0,
                        br.width(),
                        br.height(),
                    )
                )
        return chain_rects

    def fit_scene(self):
        if not self._closed:
            rect = self.scene.itemsBoundingRect().adjusted(-35, -35, 35, 35)
            self.view.fitInView(rect, QtCore.Qt.KeepAspectRatio)

    def select_nucleotide(self, index):
        self.select_indices([index], replace=True)

    def copy_dbn(self):
        text = ">%s\n%s\n%s" % (
            self.model.title,
            self.model.raw_sequence or self.model.sequence,
            self.model.raw_structure or self.model.structure,
        )
        try:
            QtWidgets.QApplication.clipboard().setText(text)
            self.status_label.setText("DBN copied to clipboard")
        except Exception as e:
            self.status_label.setText("Clipboard error: %s" % str(e))

    def export_image(self):
        filters = "PNG image (*.png)"
        if QtSvg is not None:
            filters += ";;SVG image (*.svg)"
        try:
            result = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Export RNA 2D image",
                "rna_2d.png",
                filters,
            )
            path = result[0] if isinstance(result, (tuple, list)) else result
        except Exception:
            path = ""
        if not path:
            return
        path = str(path)
        try:
            if path.lower().endswith(".svg") and QtSvg is not None:
                self._export_svg(path)
            else:
                if not path.lower().endswith(".png"):
                    path += ".png"
                self._export_png(path)
            self.status_label.setText("Exported: %s" % path)
        except Exception as e:
            msg = "Export failed: %s" % str(e)
            self.status_label.setText(msg)
            try:
                QtWidgets.QMessageBox.critical(self, "RNA 2D export", msg)
            except Exception:
                pass

    def _export_png(self, path):
        rect = self.scene.itemsBoundingRect().adjusted(-30, -30, 30, 30)
        scale = 2.0
        width = max(1, int(math.ceil(rect.width() * scale)))
        height = max(1, int(math.ceil(rect.height() * scale)))
        max_dim = 12000
        if width > max_dim or height > max_dim:
            factor = min(float(max_dim) / width, float(max_dim) / height)
            width = max(1, int(width * factor))
            height = max(1, int(height * factor))
        image = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
        image.fill(
            QtGui.QColor("#0f172a" if getattr(self, "is_dark", False) else "#ffffff")
        )
        painter = QtGui.QPainter(image)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        self.scene.render(
            painter,
            QtCore.QRectF(0, 0, width, height),
            rect,
            QtCore.Qt.KeepAspectRatio,
        )
        painter.end()
        if not image.save(path):
            raise CmdException("Qt could not save PNG file")

    def _export_svg(self, path):
        """Export the 2D RNA diagram to an SVG vector file with robust coordinate handling."""
        # 1. Obtain scene bounds with generous margin
        bounding_rect_f = self.scene.itemsBoundingRect().adjusted(-30, -30, 30, 30)

        # 2. Convert to integer QRect for cross-platform QtSvg compatibility
        view_box = bounding_rect_f.toRect()

        # Ensure non-zero width and height
        if view_box.width() <= 0:
            view_box.setWidth(100)
        if view_box.height() <= 0:
            view_box.setHeight(100)

        # 3. Configure the SVG generator
        generator = QtSvg.QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(view_box.size())
        generator.setViewBox(view_box)
        generator.setTitle(self.model.title)
        generator.setDescription("RNA secondary structure derived by DSSR")

        # 4. Render scene to SVG
        painter = QtGui.QPainter(generator)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        painter.setRenderHint(QtGui.QPainter.TextAntialiasing, True)
        try:
            self.scene.render(
                painter,
                QtCore.QRectF(view_box),
                bounding_rect_f,
                QtCore.Qt.KeepAspectRatio,
            )
        finally:
            painter.end()

    def _capture_positions(self):
        return [(float(node.pos().x()), float(node.pos().y())) for node in self.nodes]

    def _apply_positions(self, positions):
        if len(positions) != len(self.nodes):
            raise CmdException(
                "Layout has %d coordinates but this RNA has %d nucleotides"
                % (len(positions), len(self.nodes))
            )
        for node, point in zip(self.nodes, positions):
            node.setPos(float(point[0]), float(point[1]))
        self._update_scene_rect()

    def _apply_sparse_positions(self, sparse_coords):
        """Apply coordinates only to the specific nucleotide nodes that changed."""
        total = len(self.nodes)
        for index, point in sparse_coords.items():
            if 0 <= index < total:
                self.nodes[index].setPos(float(point[0]), float(point[1]))
        self._update_scene_rect()

    def _push_history(self, before, after, label="edit"):
        """Store a memory-efficient sparse diff of only the bases that changed position."""
        if not before or not after or len(before) != len(after):
            return

        diff_before = {}
        diff_after = {}

        for index, (b_pt, a_pt) in enumerate(zip(before, after)):
            if abs(b_pt[0] - a_pt[0]) > 1.0e-5 or abs(b_pt[1] - a_pt[1]) > 1.0e-5:
                diff_before[index] = (float(b_pt[0]), float(b_pt[1]))
                diff_after[index] = (float(a_pt[0]), float(a_pt[1]))

        # Nothing changed; don't add redundant undo states
        if not diff_before:
            return

        self._undo.append(
            {
                "before": diff_before,
                "after": diff_after,
                "label": str(label),
            }
        )
        if len(self._undo) > self.HISTORY_LIMIT:
            self._undo = self._undo[-self.HISTORY_LIMIT :]
        self._redo = []
        self._update_history_buttons()

    def undo_layout(self):
        if not self._undo:
            return
        entry = self._undo.pop()
        self._apply_sparse_positions(entry["before"])
        self._redo.append(entry)
        self._update_history_buttons()
        self._update_editor_status("undo: %s" % entry["label"])

    def redo_layout(self):
        if not self._redo:
            return
        entry = self._redo.pop()
        self._apply_sparse_positions(entry["after"])
        self._undo.append(entry)
        self._update_history_buttons()
        self._update_editor_status("redo: %s" % entry["label"])

    def reset_layout(self):
        before = self._capture_positions()
        algorithm = self.layout_combo.currentText().strip().lower() or "standard"
        after = Dssr2DLayout.compute(self.model, algorithm)
        self._auto_positions = [(float(x), float(y)) for x, y in after]
        self._apply_positions(after)
        self._push_history(before, self._capture_positions(), "reset layout")
        self.fit_scene()
        self._update_editor_status("automatic layout reset")

    def reset_selected_bases(self):
        if len(self._auto_positions) != len(self.nodes):
            self._auto_positions = Dssr2DLayout.compute(
                self.model,
                self.layout_combo.currentText().strip().lower() or "standard",
            )
        selected = [node for node in self.nodes if node.isSelected()]
        if not selected:
            return
        before = self._capture_positions()
        for node in selected:
            x, y = self._auto_positions[node.nt_index]
            node.setPos(float(x), float(y))
        after = self._capture_positions()
        self._push_history(before, after, "reset selected")
        self._update_scene_rect()
        self._update_editor_status("selected bases reset")

    def nudge_selected(self, dx, dy):
        selected = [node for node in self.nodes if node.isSelected()]
        if not selected:
            return
        before = self._capture_positions()
        for node in selected:
            node.moveBy(float(dx), float(dy))
        after = self._capture_positions()
        self._push_history(before, after, "nudge")
        self._update_editor_status("nudged %d base(s)" % len(selected))

    def select_indices(self, indices, replace=True):
        wanted = set(int(index) for index in indices)
        blocked = self.scene.blockSignals(True)
        try:
            for node in self.nodes:
                node.setSelected(
                    node.nt_index in wanted or (not replace and node.isSelected())
                )
        finally:
            self.scene.blockSignals(blocked)
        self._update_editor_status("selection changed")
        self._sync_pymol_selection()

    def select_all_bases(self):
        self.select_indices(range(len(self.nodes)), replace=True)

    def clear_base_selection(self):
        self.select_indices(())

    def _pair_partner(self, index):
        table = self._pair_table
        if 0 <= index < len(table):
            return int(table[index])
        return -1

    def _stem_indices(self, index):
        for stem in self._stems:
            indices = {int(i) for pair in stem.get("pairs", []) for i in pair}
            if index in indices:
                return sorted(indices)
        return self._pair_group(index)

    def _selection_changed(self):
        if self._closed or self._rebuilding:
            return
        self._schedule_live_sync()
        self._update_editor_status("selection changed")

    def _sync_pymol_selection(self):
        if self._closed:
            return
        self._sync_timer.stop()
        self._sync_pending = False
        signature = self._node_residue_signature()
        self._select_residues_in_pymol(signature)
        self._last_pymol_signature = signature
        self._update_pymol_highlight(signature)

    def save_layout(self):
        default_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.model.title)
        default_name = (default_name.strip("_") or "rna_2d") + ".dssr2d.json"
        path, _chosen = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save RNA 2D layout",
            default_name,
            "DSSR RNA 2D layout (*.dssr2d.json *.json);;All files (*)",
        )
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".dssr2d.json"
        payload = {
            "format": "DSSR-PyMOL-RNA2D",
            "version": 1,
            "title": self.model.title,
            "sequence": self.model.sequence,
            "structure": self.model.structure,
            "layout": self.layout_combo.currentText().strip().lower(),
            "positions": [
                [round(x, 6), round(y, 6)] for x, y in self._capture_positions()
            ],
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
        self._update_editor_status("layout saved: %s" % path)

    def load_layout(self):
        path, _chosen = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Load RNA 2D layout",
            "",
            "DSSR RNA 2D layout (*.dssr2d.json *.json);;All files (*)",
        )
        if not path:
            return
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        positions = payload.get("positions", [])
        if len(positions) != len(self.nodes):
            raise CmdException(
                "Saved layout contains %d bases; current RNA contains %d"
                % (len(positions), len(self.nodes))
            )
        saved_sequence = str(payload.get("sequence", ""))
        if saved_sequence and saved_sequence != self.model.sequence:
            answer = QtWidgets.QMessageBox.question(
                self,
                "Sequence mismatch",
                "This layout was saved for a different sequence. Load its "
                "coordinates anyway?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No,
            )
            if answer != QtWidgets.QMessageBox.Yes:
                return
        before = self._capture_positions()
        self._apply_positions(positions)
        self._push_history(before, self._capture_positions(), "load layout")
        self.fit_scene()
        self._update_editor_status("layout loaded: %s" % path)

    def _update_history_buttons(self):
        for button, history, action in (
            (self.undo_btn, self._undo, "Undo"),
            (self.redo_btn, self._redo, "Redo"),
        ):
            button.setEnabled(bool(history))
            if history:
                button.setToolTip("%s: %s" % (action, history[-1].get("label", "edit")))

    def _update_editor_status(self, action=""):
        selected = self._sync_sequence_selection()
        variant = str(getattr(self.model, "_dssr2d_layout_variant", self.algorithm))
        text = "Selected %d · %s · %s" % (
            len(selected),
            self.interaction_tool(),
            variant,
        )
        self.status_label.setText(text + (" · " + action if action else ""))

    def _schedule_scene_rect(self):
        if self._closed or self._scene_rect_pending:
            return
        self._scene_rect_pending = True

        def update():
            self._scene_rect_pending = False
            self._update_scene_rect()

        QtCore.QTimer.singleShot(0, update)

    def _update_scene_rect(self):
        if not self._closed:
            self.scene.setSceneRect(
                self.scene.itemsBoundingRect().adjusted(-120, -120, 120, 120)
            )

    def drag_mode(self):
        return str(self.drag_mode_combo.currentData() or "base")

    def set_drag_mode(self, mode):
        index = self.drag_mode_combo.findData(str(mode or "base").strip().lower())
        if index >= 0:
            self.drag_mode_combo.setCurrentIndex(index)

    def _pair_group(self, index):
        partner = self._pair_partner(index)
        return [index] if partner < 0 else sorted({index, partner})

    def _loop_group(self, index):
        n = len(self.model.nts)
        if index < 0 or index >= n:
            return []
        table = self._pair_table
        if table[index] >= 0:
            return self._pair_group(index)

        left = index
        right = index
        while (
            left > 0
            and (left - 1) not in self.model.chain_breaks
            and table[left - 1] < 0
        ):
            left -= 1
        while (
            right + 1 < n
            and right not in self.model.chain_breaks
            and table[right + 1] < 0
        ):
            right += 1
        result = set(range(left, right + 1))

        if left > 0 and right + 1 < n:
            a = left - 1
            b = right + 1
            if table[a] == b:
                result.add(a)
                result.add(b)
        return sorted(result)

    def _chain_segment(self, index):
        n = len(self.model.nts)
        left = 0
        right = n - 1
        for break_after in sorted(self.model.chain_breaks):
            if break_after < index:
                left = break_after + 1
            elif break_after >= index:
                right = break_after
                break
        return list(range(max(0, left), min(n - 1, right) + 1))

    def _branch_group(self, index):
        candidates = [
            stem
            for stem in self._stems
            if int(stem.get("outer_i", -1)) <= index <= int(stem.get("outer_j", -1))
        ]
        if not candidates:
            return self._chain_segment(index)
        stem = min(
            candidates,
            key=lambda item: int(item.get("outer_j", 0)) - int(item.get("outer_i", 0)),
        )
        return list(range(int(stem["outer_i"]), int(stem["outer_j"]) + 1))

    def _build_adjacency(self):
        n = len(self.model.nts)
        adjacency = [set() for _ in range(n)]
        for index in range(n - 1):
            if index not in self.model.chain_breaks:
                adjacency[index].add(index + 1)
                adjacency[index + 1].add(index)
        for pair in self.model.secondary_pairs:
            i = int(pair.get("i", -1))
            j = int(pair.get("j", -1))
            if 0 <= i < n and 0 <= j < n and i != j:
                adjacency[i].add(j)
                adjacency[j].add(i)
        self._adjacency = adjacency
        return adjacency

    def _soft_weights(self, anchors):
        anchors = sorted({int(i) for i in anchors if 0 <= int(i) < len(self.model.nts)})
        if not anchors:
            return {}
        adjacency = self._adjacency or self._build_adjacency()
        strength = float(self.follow_spin.value())
        max_depth = 3
        distance = {index: 0 for index in anchors}
        queue = deque(anchors)
        while queue:
            current = queue.popleft()
            depth = distance[current]
            if depth >= max_depth:
                continue
            for neighbor in adjacency[current]:
                if neighbor not in distance:
                    distance[neighbor] = depth + 1
                    queue.append(neighbor)
        return {
            index: 1.0 if depth == 0 else max(0.08, strength**depth)
            for index, depth in distance.items()
        }

    def _prepare_node_drag(self, node, modifiers):
        index = node.nt_index
        selected = sorted(item.nt_index for item in self.nodes if item.isSelected())
        mode = "base" if modifiers & QtCore.Qt.AltModifier else self.drag_mode()
        additive = bool(
            modifiers & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
        )
        grouped = {
            "pair": self._pair_group,
            "loop": self._loop_group,
            "stem": self._stem_indices,
            "branch": self._branch_group,
        }
        weights = None
        if additive:
            group = selected or [index]
        elif mode in grouped:
            group = grouped[mode](index)
            self.select_indices(group, replace=True)
        elif mode == "selection":
            group = selected or [index]
            if not selected:
                self.select_indices(group, replace=True)
        else:
            group = selected if index in selected and len(selected) > 1 else [index]
            if group == [index]:
                self.select_indices(group, replace=True)
            if mode == "soft":
                weights = self._soft_weights(group)
                group = sorted(weights)
        if weights is None:
            weights = dict.fromkeys(group, 1.0)
        node._drag_starts = {
            item: QtCore.QPointF(self.nodes[item].pos())
            for item in group
            if 0 <= item < len(self.nodes)
        }
        node._drag_weights = {
            item: float(weights.get(item, 1.0)) for item in node._drag_starts
        }
        self._update_editor_status("dragging %s" % mode)

    def fit_selected(self):
        selected = [node for node in self.nodes if node.isSelected()]
        if not selected:
            self.fit_scene()
            return
        rect = QtCore.QRectF(selected[0].sceneBoundingRect())
        for node in selected[1:]:
            rect = rect.united(node.sceneBoundingRect())
        self.view.fitInView(rect.adjusted(-80, -80, 80, 80), QtCore.Qt.KeepAspectRatio)

    def _ensure_animation(self):
        if not self._closed and self._view_active and not self._timer.isActive():
            self._timer.start()

    def _animation_tick(self):
        active = [node._advance_visual() for node in self.nodes]
        if not any(active):
            self._timer.stop()

    def interaction_tool(self):
        return str(self.interaction_combo.currentData() or "edit")

    def set_interaction_tool(self, tool):
        index = self.interaction_combo.findData(str(tool or "edit").lower())
        if index >= 0:
            self.interaction_combo.setCurrentIndex(index)

    def _interaction_changed(self, *_args):
        self.view.cancel_selection_gesture()
        self._update_view_cursor()
        self._update_editor_status("tool=%s" % self.interaction_tool())

    def _update_view_cursor(self):
        cursor = (
            QtCore.Qt.CrossCursor
            if self.interaction_tool() == "brush"
            else QtCore.Qt.ArrowCursor
        )
        self.view.setCursor(cursor)

    def gel_style_enabled(self):
        return self.gel_style_cb.isChecked()

    def _gel_mode_toggled(self, checked):
        self._refresh_scene_style()
        self._ensure_animation()
        self._update_editor_status("gel mode on" if checked else "gel mode off")

    def _circles_toggled(self, checked):
        self.show_circles = bool(checked)
        self._refresh_scene_style()

    def _refresh_scene_style(self):
        is_dark = getattr(self, "is_dark", False)
        label_color = (
            QtGui.QColor(248, 250, 252) if is_dark else QtGui.QColor(15, 23, 42)
        )
        for node in self.nodes:
            node.base_text_item.setBrush(
                QtGui.QBrush(
                    DssrUtils.base_text_color(node.nt.get("base", ""), self.base_colors)
                )
            )
            for child in node.childItems():
                if child is not node.base_text_item and isinstance(
                    child, QtWidgets.QGraphicsSimpleTextItem
                ):
                    child.setBrush(QtGui.QBrush(label_color))
            node.update()
        self.sequence_view.set_letter_colors(self.base_colors)
        self.scene.update()

    def _after_brush_selection(self, final=False):
        if final:
            self._flush_live_sync(final=True)
            self._update_editor_status("brush selection mapped to 3D")
        else:
            self._schedule_live_sync()
            self._update_editor_status("brushing")

    def _schedule_live_sync(self):
        if not self._closed and self._view_active:
            self._sync_pending = True
            self._sync_timer.start()

    def _flush_live_sync(self, final=False):
        if getattr(self, "_closed", False):
            return
        self._sync_pending = False
        self._sync_pymol_selection()
        if final:
            try:
                if self.zoom_3d_cb.isChecked() and cmd.count_atoms(
                    self.HIGHLIGHT_OBJECT
                ):
                    cmd.zoom(self.HIGHLIGHT_OBJECT, buffer=4.0)
            except Exception:
                pass

    def _node_residue_signature(self):
        return tuple(
            sorted(
                {
                    (
                        str(node.nt.get("chain", "")).strip(),
                        str(node.nt.get("resi", "")).strip(),
                    )
                    for node in self.nodes
                    if node.isSelected() and str(node.nt.get("resi", "")).strip()
                }
            )
        )

    def _reverse_sync_toggled(self, checked):
        self._last_pymol_signature = None
        self._last_sel_count = -1
        self._last_active_names = ()
        if checked:
            self._pull_pymol_selection()
            self._update_editor_status("bidirectional sync on")
        else:
            self._update_editor_status("3D-to-2D sync off")

    def _pull_pymol_selection(self):
        """Synchronize 3D PyMOL selections to 2D nodes without unnecessary overhead."""
        if (
            self._closed
            or self._sync_pending
            or self._sync_from_pymol
            or not self.isVisible()
            or not self.reverse_3d_cb.isChecked()
        ):
            return

        # 1. Quick check: retrieve only currently enabled selection names
        try:
            enabled_selections = tuple(cmd.get_names("selections", enabled_only=1))
        except Exception:
            enabled_selections = ()

        # 2. Fast exit if no selections exist
        if not enabled_selections:
            self._last_active_names = ()
            self._last_sel_count = 0
            if (
                any(node.isSelected() for node in self.nodes)
                and self._last_pymol_signature
            ):
                self._last_pymol_signature = tuple()
                self._sync_from_pymol = self._rebuilding = True
                try:
                    for node in self.nodes:
                        node.setSelected(False)
                finally:
                    self._rebuilding = self._sync_from_pymol = False
                self._update_pymol_highlight()
                self._update_editor_status("3D selection cleared")
            return

        # 3. Fast filter: check atom count before doing full iteration
        try:
            active_sel = " or ".join("(%s)" % s for s in enabled_selections)
            scoped = "((%s) and (%s))" % (self.pymol_selection, active_sel)
            current_count = int(cmd.count_atoms(scoped))
        except Exception:
            current_count = 0

        # If selection names and atom count have not changed, skip iteration
        if (
            enabled_selections == self._last_active_names
            and current_count == self._last_sel_count
        ):
            return

        self._last_active_names = enabled_selections
        self._last_sel_count = current_count

        if current_count <= 0:
            if self._last_pymol_signature:
                self._last_pymol_signature = tuple()
                self._sync_from_pymol = self._rebuilding = True
                try:
                    for node in self.nodes:
                        node.setSelected(False)
                finally:
                    self._rebuilding = self._sync_from_pymol = False
                self._update_pymol_highlight()
            return

        # 4. Only iterate when an actual selection change is verified
        residues = set()
        try:
            cmd.iterate(
                scoped,
                "_dssr_residues.add((chain, resi))",
                space={"_dssr_residues": residues},
            )
        except Exception:
            residues = set()

        signature = tuple(sorted(residues))
        if signature == self._last_pymol_signature:
            return
        self._last_pymol_signature = signature

        # Map residues back to 2D graph nodes
        wanted = {
            index
            for index, nt in enumerate(self.model.nts)
            if str(nt.get("resi", "")).strip()
            and (str(nt.get("chain", "")).strip(), str(nt.get("resi", "")).strip())
            in residues
        }
        if wanted == {node.nt_index for node in self.nodes if node.isSelected()}:
            return

        self._sync_from_pymol = self._rebuilding = True
        try:
            for node in self.nodes:
                node.setSelected(node.nt_index in wanted)
        finally:
            self._rebuilding = self._sync_from_pymol = False

        self._update_pymol_highlight()
        self._update_editor_status("3D selection mirrored to 2D")

    def _update_pymol_highlight(self, signature=None):
        name = self.HIGHLIGHT_OBJECT
        try:
            cmd.delete(name)
            _DSSR_BLOCK_OBJECTS.discard(name)
        except Exception:
            pass
        self._last_highlight_signature = None

    def _select_residues_in_pymol(self, signature=None):
        if signature is None:
            signature = self._node_residue_signature()
        if not signature and not any(node.isSelected() for node in self.nodes):
            try:
                cmd.select("sele", "none")
            except Exception:
                pass
            return

        core = DssrParser._compact_sel_from_residues(
            {(chain, resi) for chain, resi in signature if chain}
        )
        pieces = ["(%s)" % core] if core else []
        pieces.extend("(resi %s)" % resi for chain, resi in signature if not chain)
        if pieces:
            expression = "((%s) and (%s))" % (
                self.pymol_selection,
                " or ".join(pieces),
            )
            try:
                cmd.select("sele", "byres (%s)" % expression)
            except Exception as error:
                self.status_label.setText("PyMOL selection error: %s" % str(error))

    def _build_widgets(self):
        """Build the toolbars, sequence view, and single 2D graphics canvas."""
        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        top = QtWidgets.QHBoxLayout()
        root.addLayout(top)
        self.layout_combo = DssrUI.combo(LAYOUT_CHOICES)
        self.layout_combo.setCurrentText(self.algorithm)
        top.addWidget(self.layout_combo)
        self.fit_btn = DssrUI.button("Fit", self.fit_scene, "Fit drawing [F]")
        self.undo_btn = DssrUI.button("Undo", self.undo_layout, "Undo [Ctrl+Z]")
        self.redo_btn = DssrUI.button("Redo", self.redo_layout, "Redo [Ctrl+Y]")
        for button in (self.fit_btn, self.undo_btn, self.redo_btn):
            top.addWidget(button)
        self.redraw_btn = DssrUI.button("Reset layout", self.reset_layout)
        top.addWidget(self.redraw_btn)
        top.addStretch(1)
        self.file_menu_btn = DssrUI.menu_button(
            "File",
            (
                ("Export image…", self.export_image),
                ("Save layout…", self.save_layout),
                ("Load layout…", self.load_layout),
                ("Copy sequence / DBN", self.copy_dbn),
            ),
            self,
        )
        top.addWidget(self.file_menu_btn)
        self.selection_menu_btn = DssrUI.menu_button(
            "Selection",
            (
                ("Select all [Ctrl+A]", self.select_all_bases),
                ("Clear selection [Esc]", self.clear_base_selection),
                ("Fit selected [C]", self.fit_selected),
                ("Reset selected coordinates", self.reset_selected_bases),
            ),
            self,
        )
        top.addWidget(self.selection_menu_btn)
        self.options_btn = QtWidgets.QPushButton("Options")
        self.options_btn.setCheckable(True)
        top.addWidget(self.options_btn)

        tools = QtWidgets.QHBoxLayout()
        root.addLayout(tools)
        self.interaction_combo = QtWidgets.QComboBox()
        self.interaction_combo.addItem("Select / edit [P]", "edit")
        self.interaction_combo.addItem("Brush select [B]", "brush")
        self.interaction_combo.currentIndexChanged.connect(self._interaction_changed)
        tools.addWidget(self.interaction_combo)
        self.brush_radius_spin = DssrUI.spinbox(
            8, 100, 32, suffix=" px", tip="Brush radius"
        )
        tools.addWidget(self.brush_radius_spin)
        self.drag_mode_combo = QtWidgets.QComboBox()
        for text, value in (
            ("Soft drag", "soft"),
            ("Single base", "base"),
            ("Selected bases", "selection"),
            ("Base pair", "pair"),
            ("Loop", "loop"),
            ("Stem", "stem"),
            ("Branch", "branch"),
        ):
            self.drag_mode_combo.addItem(text, value)
        self.drag_mode_combo.currentIndexChanged.connect(
            lambda: self._update_editor_status("drag mode changed")
        )
        tools.addWidget(self.drag_mode_combo)
        self.gel_style_cb = DssrUI.checkbox(
            "Gel",
            False,
            self._gel_mode_toggled,
            "Glass-like rendering and elastic motion",
        )
        tools.addWidget(self.gel_style_cb)
        tools.addStretch(1)

        self.options_panel = QtWidgets.QGroupBox("Display and 3D")
        options = QtWidgets.QGridLayout(self.options_panel)
        self.number_spin = DssrUI.spinbox(0, 10000, self.number_every)
        self.tertiary_cb = DssrUI.checkbox("Extra DSSR pairs", self.show_tertiary)
        self.base_colors_cb = DssrUI.checkbox(
            "Base colors", True, tip="Color bases using classical PyMOL/DSSR scheme"
        )
        self.circles_cb = DssrUI.checkbox(
            "Circles",
            True,
            changed=self._circles_toggled,
            tip="Show or hide circular node borders around bases",
        )
        self.follow_spin = DssrUI.spinbox(0.1, 0.9, 0.62, decimals=True, step=0.05)
        self.live_3d_cb = DssrUI.checkbox(
            "3D highlight", False, self._sync_pymol_selection
        )
        self.reverse_3d_cb = DssrUI.checkbox(
            "3D → 2D sync", True, self._reverse_sync_toggled
        )
        self.zoom_3d_cb = DssrUI.checkbox("Zoom after brush", False)
        options.addWidget(QtWidgets.QLabel("Number every"), 0, 0)
        options.addWidget(self.number_spin, 0, 1)
        options.addWidget(self.tertiary_cb, 0, 2)
        options.addWidget(self.base_colors_cb, 0, 3)
        options.addWidget(self.circles_cb, 0, 4)
        options.addWidget(QtWidgets.QLabel("Elasticity"), 1, 0)
        options.addWidget(self.follow_spin, 1, 1)
        options.addWidget(self.live_3d_cb, 1, 2)
        options.addWidget(self.reverse_3d_cb, 1, 3)
        options.addWidget(self.zoom_3d_cb, 2, 2, 1, 2)
        root.addWidget(self.options_panel)
        self.options_panel.hide()
        self.options_btn.toggled.connect(self.options_panel.setVisible)

        self.sequence_view = Dssr2DSequenceView(self)
        root.addWidget(self.sequence_view)
        self.scene = QtWidgets.QGraphicsScene(self)
        self.view = Dssr2DGraphicsView(self.scene, self)
        root.addWidget(self.view, 1)
        hint = QtWidgets.QLabel(
            "Drag blank space: box select · Drag bases: edit · Arrows / WASD: pan · Wheel: zoom · B: brush"
        )
        hint.setObjectName("studioHint")
        hint.setWordWrap(True)
        root.addWidget(hint)
        self.status_label = QtWidgets.QLabel(self.model.summary())
        self.status_label.setWordWrap(True)
        root.addWidget(self.status_label)

    def set_view_active(self, active):
        active = bool(active) and not self._closed
        if active == self._view_active:
            return
        self._view_active = active
        if active:
            self._last_pymol_signature = None
            self._pull_pymol_selection()
            self._reverse_timer.start()
            self._ensure_animation()
        else:
            self.view.cancel_selection_gesture()
            self.sequence_view._drag_anchor = None
            self._sync_pending = False
            for timer in (self._timer, self._sync_timer, self._reverse_timer):
                timer.stop()
            cmd.delete(self.HIGHLIGHT_OBJECT)
            _DSSR_BLOCK_OBJECTS.discard(self.HIGHLIGHT_OBJECT)
            self._last_highlight_signature = None

    def showEvent(self, event):
        super().showEvent(event)
        if not self._shown_once:
            self._shown_once = True
            QtCore.QTimer.singleShot(0, self.fit_scene)

    def shutdown(self):
        if not self._closed:
            self.set_view_active(False)
            self._closed = True
            self.setEnabled(False)

    def closeEvent(self, event):
        self.shutdown()
        super().closeEvent(event)

    def _sync_sequence_selection(self):
        selected = [node.nt_index for node in self.nodes if node.isSelected()]
        self.sequence_view.set_selected_indices(selected)
        return selected


# Public PyMOL commands are registered once
dssr_select = DssrCmd.dssr_select
dssr_block = DssrCmd.dssr_block
dssr_2d = DssrCmd.dssr_2d
dssr_gui = DssrGuiDialog.dssr_gui

for _command in (dssr_select, dssr_block, dssr_2d, dssr_gui):
    cmd.extend(_command.__name__, _command)

try:
    for _name in ("dssr_select", "dssr_block", "dssr_2d"):
        cmd.auto_arg[0][_name] = cmd.auto_arg[0]["zoom"]
    cmd.auto_arg[2]["dssr_block"] = [cmd.Shortcut(BLOCK_FEATURES), "block_file", ""]
except (AttributeError, KeyError):
    pass


def __init_plugin__(app=None):
    addmenuitemqt("DSSR", dssr_gui)


print("Loaded DSSR-PyMOL %s (RNA 2D studio, Python/Qt)" % __DSSR_PLUGIN_VERSION__)
