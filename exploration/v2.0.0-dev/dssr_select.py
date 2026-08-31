# dssr_select.py
# DSSR-based selection of RNA structural features in PyMOL
#
# (c) 2026 Bener Dulger, Eric Chen, and Xiang-Jun Lu
#
# This project was initiated and coordinated by Xiang-Jun Lu.
#
# CONTRIBUTIONS:
# - Bener Dulger: Initial implementation of structural feature selection and JSON parsing;
#                 engineered the modular, class-based architectural refactor for v1.1.0-dev.
# - Eric Chen: Developed the Qt-based GUI, incorporated dssr_block functionality,
#              enhanced feature parsing, and performed code consolidation of v1.0.0.
# - Thomas Holder: Original 'dssr_block' logic (c) Schrodinger LLC.
#
# LICENSE: BSD 2-Clause
#
# This plugin incorporates code from 'dssr_block'.
# Redistributions must retain the original copyright notice and this license.

from pymol import cmd, CmdException
from pymol.Qt import QtWidgets, QtCore
from pymol.plugins import addmenuitemqt
import subprocess
import json
import re
import tempfile
import os

__DSSR_PLUGIN_VERSION__ = "v1.3.3-white-ui-2026-08-18"
_DSSR_GUI_DIALOG = None
_hex_color_cache = {}
selected_features = []
_DSSR_BLOCK_OBJECTS = set()
_DSSR_SELECTION_OBJECTS = set()

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


class HelperFunctions:
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
    def run_dssr_json(pdb_path, exe):

        # Keep DSSR identifiers compatible with Jmol/EBI unit IDs and request
        # the U-turn annotations exposed by the current upstream plugin.
        args = [exe, "--json", "--u-turn", "--idstr=ebi", "-i=" + pdb_path]

        try:
            p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            rc = p.returncode
        except OSError:
            raise CmdException('Cannot execute exe="%s"' % exe)

        try:
            out_txt = out.decode("utf-8", errors="replace") if out else ""
        except Exception:
            out_txt = str(out)
        try:
            err_txt = err.decode("utf-8", errors="replace") if err else ""
        except Exception:
            err_txt = str(err)

        if rc != 0:
            raise CmdException(
                "DSSR failed (rc=%s). stderr tail: %s"
                % (str(rc), HelperFunctions._safe_tail(err_txt))
            )

        if not out_txt.strip():
            raise CmdException(
                "DSSR returned empty stdout (expected JSON). stderr tail: %s"
                % HelperFunctions._safe_tail(err_txt)
            )

        try:
            return json.loads(out_txt)
        except Exception:
            s = out_txt
            i = s.find("{")
            j = s.rfind("}")
            if i >= 0 and j > i:
                try:
                    return json.loads(s[i : j + 1])
                except Exception as e2:
                    raise CmdException(
                        "Failed to parse DSSR JSON. stdout head: %s | stderr tail: %s | err: %s"
                        % (
                            s[:120].replace("\n", " "),
                            HelperFunctions._safe_tail(err_txt),
                            str(e2),
                        )
                    )
            raise CmdException(
                "Failed to parse DSSR JSON (no JSON object found). stdout head: %s | stderr tail: %s"
                % (s[:120].replace("\n", " "), HelperFunctions._safe_tail(err_txt))
            )

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
            hex6, rgb = HelperFunctions._hex_to_rgb01(s)
            if hex6 in _hex_color_cache:
                return _hex_color_cache[hex6]
            cname = "dssr_hex_%s" % hex6
            cmd.set_color(cname, rgb)
            _hex_color_cache[hex6] = cname
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


class ParsingAlgos:
    @staticmethod
    def feature_entries(dssr_data, feature):
        """Return a normalized list for a FEATURE_MAP entry.

        DSSR normally emits arrays, but ``nonStack`` can be emitted as a
        single object.  Keeping that compatibility rule here prevents the GUI,
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
        """
        Parse pseudoknot layers from DSSR/Vienna dot-bracket notation.

        Canonical ``()`` pairs are intentionally excluded.  Square/curly/
        angle brackets and extended ``A...a`` notation are returned as
        pseudoknot layers.  Chain separators do not consume a nucleotide
        index, keeping the resulting indices aligned with DSSR's ``nts``
        array for multi-chain structures.
        """
        stack_map = {"(": ")", "[": "]", "{": "}", "<": ">"}
        close_to_open = {v: k for k, v in stack_map.items()}
        stacks = {}
        layers = {}
        layer_assignment = {"()": 0, "[]": 1, "{}": 2, "<>": 3}
        letter_layers = {}
        next_layer = 4

        nt_index = 0
        for char in str(dotbracket):
            if char.isspace() or char in ("&", "+"):
                continue

            if char in stack_map:
                bracket_type = char + stack_map[char]
                stacks.setdefault(bracket_type, []).append(nt_index)

            elif char in close_to_open:
                open_char = close_to_open[char]
                bracket_type = open_char + char
                stack = stacks.setdefault(bracket_type, [])
                if stack:
                    open_idx = stack.pop()
                    layer = layer_assignment[bracket_type]
                    if layer != 0:
                        layers.setdefault(layer, []).append((open_idx, nt_index))

            elif "A" <= char <= "Z":
                stacks.setdefault(char, []).append(nt_index)
                if char not in letter_layers:
                    letter_layers[char] = next_layer
                    next_layer += 1

            elif "a" <= char <= "z":
                opener = char.upper()
                stack = stacks.setdefault(opener, [])
                if stack:
                    open_idx = stack.pop()
                    layer = letter_layers.get(opener)
                    if layer is None:
                        layer = next_layer
                        letter_layers[opener] = layer
                        next_layer += 1
                    layers.setdefault(layer, []).append((open_idx, nt_index))

            nt_index += 1

        for pairs in layers.values():
            pairs.sort()

        return layers

    @staticmethod
    def build_selection_from_layer(layer_pairs, nts_list):
        residues = set()
        for open_idx, close_idx in layer_pairs:
            if open_idx < len(nts_list):
                nt1_id = nts_list[open_idx].get("nt_id")
                if nt1_id:
                    residues.add(ParsingAlgos.parse_nt_id(nt1_id))
            if close_idx < len(nts_list):
                nt2_id = nts_list[close_idx].get("nt_id")
                if nt2_id:
                    residues.add(ParsingAlgos.parse_nt_id(nt2_id))

        if not residues:
            return None

        return " or ".join("(chain %s and resi %s)" % (c, r) for c, r in residues)

    @staticmethod
    def build_selection_from_pair(pair_entry):
        nt1 = pair_entry.get("nt1")
        nt2 = pair_entry.get("nt2")
        if not nt1 or not nt2:
            raise CmdException("Pair entry missing nt1 or nt2")
        residues = {ParsingAlgos.parse_nt_id(nt1), ParsingAlgos.parse_nt_id(nt2)}
        return " or ".join("(chain %s and resi %s)" % (c, r) for c, r in residues)

    @staticmethod
    def build_selection_from_nts_list(nts_list):
        if not nts_list:
            raise CmdException("Empty nucleotide list")
        residues = {ParsingAlgos.parse_nt_id(nt) for nt in nts_list}
        return " or ".join("(chain %s and resi %s)" % (c, r) for c, r in residues)

    @staticmethod
    def build_selection_from_stem(stem_entry):
        pairs = stem_entry.get("pairs", [])
        if not pairs:
            raise CmdException("Stem has no pairs")

        residues = set()
        for p in pairs:
            if p.get("nt1"):
                residues.add(ParsingAlgos.parse_nt_id(p["nt1"]))
            if p.get("nt2"):
                residues.add(ParsingAlgos.parse_nt_id(p["nt2"]))

        return " or ".join("(chain %s and resi %s)" % (c, r) for c, r in residues)

    @staticmethod
    def build_selection_from_hairpin(hairpin_entry):
        nts_long = hairpin_entry.get("nts_long")
        if not nts_long:
            raise CmdException("Hairpin missing nts_long field")
        nts_list = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
        return ParsingAlgos.build_selection_from_nts_list(nts_list)

    @staticmethod
    def build_selection_from_coaxstack(coax_entry, stems_list):
        stem_indices = coax_entry.get("stem_indices", [])
        if not stem_indices:
            raise CmdException("coaxStacks entry missing stem_indices")

        residues = set()
        for si in stem_indices:
            try:
                idx = int(si)
            except Exception:
                continue
            if idx < 1 or idx > len(stems_list):
                continue
            stem_entry = stems_list[idx - 1]
            pairs = stem_entry.get("pairs", [])
            for p in pairs:
                if p.get("nt1"):
                    residues.add(ParsingAlgos.parse_nt_id(p["nt1"]))
                if p.get("nt2"):
                    residues.add(ParsingAlgos.parse_nt_id(p["nt2"]))

        if not residues:
            raise CmdException("Could not build selection for coaxStacks entry")

        return " or ".join("(chain %s and resi %s)" % (c, r) for c, r in residues)

    @staticmethod
    def build_selection_from_atom2base(a2b_entry):
        atom = a2b_entry.get("atom")
        nt = a2b_entry.get("nt")

        clauses = []

        if nt:
            c_nt, r_nt = ParsingAlgos.parse_nt_id(nt)
            clauses.append("(chain %s and resi %s)" % (c_nt, r_nt))

        if atom:
            c_a, r_a, atom_name = ParsingAlgos.parse_a2b_atom(atom)
            atom_name = str(atom_name).replace('"', '\\"')
            clauses.append(
                '(chain %s and resi %s and name "%s")' % (c_a, r_a, atom_name)
            )

        if not clauses:
            raise CmdException("atom2bases entry missing atom and nt")
        return " or ".join(clauses)

    @staticmethod
    def build_selection_from_aminor(aminor_entry):
        desc_long = aminor_entry.get("desc_long", "")
        if not desc_long or "vs" not in desc_long:
            raise CmdException("Aminors entry missing desc_long")

        left, right = desc_long.split("vs", 1)
        nts = []
        left = left.strip()
        right = right.strip()

        if left:
            nts.append(left)
        if right:
            for item in right.split(","):
                item = item.strip()
                if item:
                    nts.append(item)

        if not nts:
            raise CmdException("Aminors entry has empty residues")

        return ParsingAlgos.build_selection_from_nts_list(nts)

    @staticmethod
    def build_selection_from_gquadruplex(gquad_entry):
        nts_long = gquad_entry.get("nts_long", "")
        nts_list = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
        return ParsingAlgos.build_selection_from_nts_list(nts_list)

    @staticmethod
    def build_selection_from_uturns(uturn_entry):
        nts_long = uturn_entry.get("nts_long", "")
        nts_list = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
        return ParsingAlgos.build_selection_from_nts_list(nts_list)

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
            return "%d: %s - %s%s" % (i, nt1, nt2, (" (%s)" % lw) if lw else "")

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
            s = ParsingAlgos._shorten_nts_long(entry.get("nts_long", ""))
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
        try:
            nts = dssr_data.get("nts", [])
        except Exception:
            nts = []
        if isinstance(nts, list):
            for nt in nts:
                try:
                    nt_id = nt.get("nt_id", "")
                except Exception:
                    nt_id = ""
                if nt_id:
                    try:
                        c, _ = ParsingAlgos.parse_nt_id(nt_id)
                    except Exception:
                        c = None
                    if c:
                        chains.add(str(c))
        return sorted(chains)

    @staticmethod
    def _count_pseudoknot_layers(dssr_data):
        try:
            dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
            layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket)
            return len(layers) if layers else 0
        except Exception:
            return 0

    @staticmethod
    def _format_rna_summary_text(dssr_data):
        chains = ParsingAlgos._extract_chain_names(dssr_data)
        chains_txt = " ".join(chains) if chains else "(unknown)"

        pairs_n = (
            len(dssr_data.get("pairs", []))
            if isinstance(dssr_data.get("pairs", None), list)
            else 0
        )
        hairpins_n = (
            len(dssr_data.get("hairpins", []))
            if isinstance(dssr_data.get("hairpins", None), list)
            else 0
        )
        stems_n = (
            len(dssr_data.get("stems", []))
            if isinstance(dssr_data.get("stems", None), list)
            else 0
        )
        bulges_n = (
            len(dssr_data.get("bulges", []))
            if isinstance(dssr_data.get("bulges", None), list)
            else 0
        )
        junctions_n = (
            len(dssr_data.get("junctions", []))
            if isinstance(dssr_data.get("junctions", None), list)
            else 0
        )
        pk_n = ParsingAlgos._count_pseudoknot_layers(dssr_data)
        aminors_n = (
            len(dssr_data.get("Aminors", []))
            if isinstance(dssr_data.get("Aminors", None), list)
            else 0
        )
        stacks_n = (
            len(dssr_data.get("stacks", []))
            if isinstance(dssr_data.get("stacks", None), list)
            else 0
        )

        gquads_n = (
            len(dssr_data.get("Gtetrads", []))
            if isinstance(dssr_data.get("Gtetrads", None), list)
            else 0
        )
        uturns_n = (
            len(dssr_data.get("Uturns", []))
            if isinstance(dssr_data.get("Uturns", None), list)
            else 0
        )
        lines = []
        lines.append("RNA Structure Summary")
        lines.append("---------------------")
        lines.append("Chains: %s" % chains_txt)
        lines.append("Base pairs: %d" % pairs_n)
        lines.append("Hairpins: %d" % hairpins_n)
        lines.append("Stems: %d" % stems_n)
        lines.append("Bulges: %d" % bulges_n)
        lines.append("Junctions: %d" % junctions_n)
        lines.append("Pseudoknots: %d" % pk_n)
        lines.append("A-minor interactions: %d" % aminors_n)
        lines.append("Stacking interactions: %d" % stacks_n)
        lines.append("G-quadruplexes: %d" % gquads_n)
        lines.append("U-turns: %d" % uturns_n)
        return "\n".join(lines)

    @staticmethod
    def _collect_residues_from_nts_long(nts_long):
        residues = set()
        if not nts_long:
            return residues
        nts_list = [nt.strip() for nt in str(nts_long).split(",") if nt.strip()]
        for nt in nts_list:
            try:
                residues.add(ParsingAlgos.parse_nt_id(nt))
            except Exception:
                pass
        return residues

    @staticmethod
    def _collect_residues_all(dssr_data, feature):
        feature = str(feature).lower().strip()
        residues = set()

        if feature == "pairs":
            json_key = FEATURE_MAP.get(feature, feature)
            pairs = dssr_data.get(json_key, [])
            if isinstance(pairs, list):
                for p in pairs:
                    try:
                        nt1 = p.get("nt1")
                        nt2 = p.get("nt2")
                    except Exception:
                        nt1, nt2 = None, None
                    if nt1:
                        try:
                            residues.add(ParsingAlgos.parse_nt_id(nt1))
                        except Exception:
                            pass
                    if nt2:
                        try:
                            residues.add(ParsingAlgos.parse_nt_id(nt2))
                        except Exception:
                            pass
            return residues

        if feature == "stems":
            stems = dssr_data.get("stems", [])
            if isinstance(stems, list):
                for st in stems:
                    try:
                        pairs = st.get("pairs", [])
                    except Exception:
                        pairs = []
                    if not isinstance(pairs, list):
                        continue
                    for p in pairs:
                        try:
                            cmd_nt1 = p.get("nt1")
                            cmd_nt2 = p.get("nt2")
                        except Exception:
                            cmd_nt1, cmd_nt2 = None, None
                        if cmd_nt1:
                            try:
                                residues.add(ParsingAlgos.parse_nt_id(cmd_nt1))
                            except Exception:
                                pass
                        if cmd_nt2:
                            try:
                                residues.add(ParsingAlgos.parse_nt_id(cmd_nt2))
                            except Exception:
                                pass
            return residues

        if feature in (
            "hairpins",
            "bulges",
            "junctions",
            "stacks",
            "nonstack",
            "gquadruplexes",
            "uturns",
        ):
            entries = ParsingAlgos.feature_entries(dssr_data, feature)
            for e in entries:
                try:
                    nts_long = e.get("nts_long", "")
                except Exception:
                    nts_long = ""
                residues |= ParsingAlgos._collect_residues_from_nts_long(nts_long)
            return residues

        if feature == "aminors":
            entries = dssr_data.get("Aminors", [])
            if isinstance(entries, list):
                for a in entries:
                    try:
                        desc_long = a.get("desc_long", "")
                    except Exception:
                        desc_long = ""
                    if not desc_long or "vs" not in desc_long:
                        continue
                    try:
                        left, right = desc_long.split("vs", 1)
                    except Exception:
                        continue
                    nts = []
                    left = left.strip()
                    right = right.strip()
                    if left:
                        nts.append(left)
                    if right:
                        for item in right.split(","):
                            item = item.strip()
                            if item:
                                nts.append(item)
                    for nt in nts:
                        try:
                            residues.add(ParsingAlgos.parse_nt_id(nt))
                        except Exception:
                            pass
            return residues

        if feature == "pseudoknot":
            try:
                dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
                nts_list = dssr_data.get("nts", None)
                if nts_list is None or not isinstance(nts_list, list):
                    return residues
                layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket) or {}
                for pairs in layers.values():
                    for open_idx, close_idx in pairs:
                        if open_idx < len(nts_list):
                            nt1_id = nts_list[open_idx].get("nt_id")
                            if nt1_id:
                                try:
                                    residues.add(ParsingAlgos.parse_nt_id(nt1_id))
                                except Exception:
                                    pass
                        if close_idx < len(nts_list):
                            nt2_id = nts_list[close_idx].get("nt_id")
                            if nt2_id:
                                try:
                                    residues.add(ParsingAlgos.parse_nt_id(nt2_id))
                                except Exception:
                                    pass
            except Exception:
                pass
            return residues

        return residues

    @staticmethod
    def _sort_resi_key(resi_str):
        try:
            return (0, int(str(resi_str)))
        except Exception:
            return (1, str(resi_str))

    @staticmethod
    def _compact_sel_from_residues(residues):
        if not residues:
            return ""
        by_chain = {}
        for c, r in residues:
            c = str(c)
            r = str(r)
            by_chain.setdefault(c, set()).add(r)

        parts = []
        for c in sorted(by_chain.keys()):
            resis = sorted(by_chain[c], key=ParsingAlgos._sort_resi_key)
            resi_expr = "+".join(resis)
            parts.append("(chain %s and resi %s)" % (c, resi_expr))
        return " or ".join(parts)

    @staticmethod
    def _build_residue_sel_from_dssr(dssr_data, feature, index):
        feature = str(feature).lower().strip()
        idx = int(index)

        if feature == "pseudoknot":
            dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
            nts_list = dssr_data.get("nts", None)
            if nts_list is None:
                raise CmdException("No nts found in DSSR output")

            layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket)
            if not layers:
                raise CmdException("No pseudoknot layers found")

            layer_keys = sorted(layers.keys())
            if idx < 1 or idx > len(layer_keys):
                raise CmdException(
                    "pseudoknot layer index %d out of range (1..%d)"
                    % (idx, len(layer_keys))
                )

            pairs = layers[layer_keys[idx - 1]]
            sel_str = ParsingAlgos.build_selection_from_layer(pairs, nts_list)
            if not sel_str:
                raise CmdException(
                    "Could not build selection for pseudoknot layer %d" % idx
                )
            return sel_str

        if feature not in FEATURE_MAP:
            raise CmdException('Unknown feature "%s"' % feature)

        json_key = FEATURE_MAP[feature]
        feature_list = ParsingAlgos.feature_entries(dssr_data, feature)
        if not feature_list:
            raise CmdException('No "%s" found in DSSR output' % json_key)

        if idx < 1 or idx > len(feature_list):
            raise CmdException(
                "%s index %d out of range (1..%d)" % (feature, idx, len(feature_list))
            )

        entry = feature_list[idx - 1]

        if feature == "pairs":
            return ParsingAlgos.build_selection_from_pair(entry)

        if feature in ("stems", "helices"):
            return ParsingAlgos.build_selection_from_stem(entry)

        if feature == "hairpins":
            return ParsingAlgos.build_selection_from_hairpin(entry)

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
        ):
            nts_long = entry.get("nts_long", "")
            if not nts_long:
                raise CmdException("%s entry missing nts_long field" % feature)
            nts_list_parsed = [nt.strip() for nt in nts_long.split(",") if nt.strip()]
            return ParsingAlgos.build_selection_from_nts_list(nts_list_parsed)

        if feature == "coaxstacks":
            stems_list = dssr_data.get("stems", [])
            if not stems_list:
                raise CmdException("No stems found, required for coaxStacks")
            return ParsingAlgos.build_selection_from_coaxstack(entry, stems_list)

        if feature == "atom2bases":
            return ParsingAlgos.build_selection_from_atom2base(entry)

        if feature == "aminors":
            return ParsingAlgos.build_selection_from_aminor(entry)

        if feature == "gquadruplexes":
            return ParsingAlgos.build_selection_from_gquadruplex(entry)

        if feature == "uturns":
            return ParsingAlgos.build_selection_from_uturns(entry)

        if feature == "nts":
            nt_id = entry.get("nt_id")
            if not nt_id:
                raise CmdException("Nucleotide entry missing nt_id field")
            c, r = ParsingAlgos.parse_nt_id(nt_id)
            return "(chain %s and resi %s)" % (c, r)

        raise CmdException('Feature "%s" not supported for residue selection' % feature)


class DssrFunctions:
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
        feature = HelperFunctions.unquote(feature).lower().strip()

        user_color = HelperFunctions._resolve_color_spec(
            HelperFunctions.unquote(color).strip()
        )

        if feature in ("features", "help"):
            keys = sorted(FEATURE_MAP.keys())
            print("Supported features: " + ", ".join(keys))
            print("Tip: use index=0 to list detected items for a feature.")
            return

        if feature not in FEATURE_MAP:
            valid = ", ".join(sorted(FEATURE_MAP.keys()))
            raise CmdException('Unknown feature "%s". Valid: %s' % (feature, valid))

        json_key = FEATURE_MAP[feature]

        if state == 0 or state < 0:
            state = cmd.get_state()

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmpfilepdb = tmp.name
        tmp.close()

        try:
            cmd.save(tmpfilepdb, selection, state)
            if precolor:
                cmd.color("gray", selection)

            dssr_data = HelperFunctions.run_dssr_json(tmpfilepdb, exe)

            if feature == "pseudoknot":
                layer_colors = ["blue", "pink", "green", "yellow", "orange"]

                dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
                nts_list = dssr_data.get("nts", None)
                if nts_list is None:
                    raise CmdException("No nts found in DSSR output")

                layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket)
                if not layers:
                    raise CmdException("No pseudoknot layers found in structure")

                layer_keys = sorted(layers.keys())

                if index == 0:
                    print("pseudoknot: %d layer(s)" % len(layer_keys))
                    for j, k in enumerate(layer_keys, 1):
                        print(
                            "  layer %d (key=%s): %d pair(s)"
                            % (j, str(k), len(layers[k]))
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

                sel_str = ParsingAlgos.build_selection_from_layer(pairs, nts_list)
                if sel_str is None:
                    raise CmdException("Could not build selection for layer %d" % index)

                scoped_sel = "((%s) and (%s))" % (selection, sel_str)
                cmd.select(name, scoped_sel)
                selected_count = int(cmd.count_atoms(name))
                if selected_count <= 0:
                    cmd.delete(name)
                    raise CmdException(
                        "DSSR residues did not map back to the requested PyMOL selection"
                    )
                _DSSR_SELECTION_OBJECTS.add(str(name))
                cmd.color(user_color if user_color else layer_color, name)

                if not quiet:
                    print(
                        'dssr_select: selection "%s" pseudoknot layer %d with %d pair(s)'
                        % (name, index, len(pairs))
                    )
                return

            feature_list = ParsingAlgos.feature_entries(dssr_data, feature)
            if not feature_list:
                raise CmdException('No "%s" found in DSSR output' % json_key)

            if index == 0:
                total = len(feature_list)
                print("%s: %d item(s)" % (feature, total))
                show_n = 20 if not quiet else 10
                show_n = min(show_n, total)
                for i in range(show_n):
                    print(
                        "  "
                        + ParsingAlgos._preview_entry(feature, feature_list[i], i + 1)
                    )
                if total > show_n:
                    print("  ... (%d more)" % (total - show_n))
                return

            if index < 1 or index > len(feature_list):
                raise CmdException(
                    "%s index %d out of range (1..%d)"
                    % (feature, index, len(feature_list))
                )

            sel_str = ParsingAlgos._build_residue_sel_from_dssr(
                dssr_data, feature, index
            )
            if not sel_str:
                raise CmdException(
                    "Could not build selection for %s index %d" % (feature, index)
                )
            scoped_sel = "((%s) and (%s))" % (selection, sel_str)
            cmd.select(name, scoped_sel)
            selected_count = int(cmd.count_atoms(name))
            if selected_count <= 0:
                cmd.delete(name)
                raise CmdException(
                    "DSSR residues did not map back to the requested PyMOL selection"
                )
            _DSSR_SELECTION_OBJECTS.add(str(name))
            cmd.color(user_color if user_color else "pink", name)
            selected_features.append(feature)

            if not quiet:
                print(
                    'dssr_select: created selection "%s" for %s (index %d) in state %d'
                    % (name, feature, index, state)
                )

        finally:
            try:
                os.remove(tmpfilepdb)
            except OSError:
                pass

    @staticmethod
    def _dssr_default_selection():
        objs = cmd.get_object_list("enabled")
        if len(objs) == 1:
            return objs[0]
        return "all"

    @staticmethod
    def dssr(
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
        selection = selection or sel or DssrFunctions._dssr_default_selection()

        feature_in = f if f is not None else feature
        feature_in = HelperFunctions.unquote(feature_in).strip()

        if feature_in.lower() in ("features", "help"):
            DssrFunctions.dssr_select(
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

        DssrFunctions.dssr_select(
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

        if int(display):
            cmd.show("sticks", nm)
            try:
                cmd.set("stick_radius", float(stick_radius), nm)
            except Exception:
                pass
            if int(do_zoom):
                cmd.zoom(nm)

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

        tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmp_r3d = tempfile.NamedTemporaryFile(suffix=".r3d", delete=False)
        tmpfilepdb = tmp_pdb.name
        tmpfiler3d = tmp_r3d.name
        tmp_pdb.close()
        tmp_r3d.close()

        if not name:
            name = DssrFunctions._unused_name("dssr_block")

        try:
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
                    "--block-file=" + HelperFunctions.unquote(block_file),
                    "--block-depth=" + str(block_depth),
                    "-i=" + tmpfilepdb,
                    "-o=" + tmpfiler3d,
                ]

                # Incorporate block_color argument
                if block_color:
                    args.append("--block-color=" + HelperFunctions.unquote(block_color))

                try:
                    p = subprocess.Popen(
                        args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    out, err = p.communicate()
                    rc = p.returncode
                except OSError:
                    raise CmdException('Cannot execute exe="%s"' % exe)

                if rc != 0:
                    err_txt = ""
                    try:
                        err_txt = err.decode("utf-8", errors="replace") if err else ""
                    except Exception:
                        err_txt = str(err)
                    raise CmdException(
                        "DSSR block failed (rc=%s). stderr tail: %s"
                        % (str(rc), HelperFunctions._safe_tail(err_txt))
                    )

                cmd.load(tmpfiler3d, name, max(1, st), zoom=0)

            _DSSR_BLOCK_OBJECTS.add(str(name))
            if not quiet:
                print(
                    'dssr_block: loaded "%s" (block_file=%s, block_depth=%s)'
                    % (name, str(block_file), str(block_depth))
                )

        finally:
            try:
                os.remove(tmpfilepdb)
            except OSError:
                pass
            try:
                os.remove(tmpfiler3d)
            except OSError:
                pass

    @staticmethod
    def _wrap_seq(s, width):
        try:
            width = int(width)
        except Exception:
            width = 80
        if width <= 0:
            return s
        return "\n".join(s[i : i + width] for i in range(0, len(s), width))

    @staticmethod
    def _revcomp(seq):
        s = "".join([c for c in str(seq).upper() if c.isalpha()])
        if "U" in s and "T" not in s:
            comp = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
        else:
            comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(comp.get(b, "N") for b in s[::-1])

    @staticmethod
    def _parse_fastastr(fasta_text):
        blocks = []
        header = None
        seq = []
        for line in str(fasta_text).splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    blocks.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            blocks.append((header, "".join(seq)))
        return blocks

    @staticmethod
    def dssr_seq(selection="all", chain="", fmt="raw", wrap=80, rc=0, quiet=1):
        fmt = str(fmt).strip().lower()
        quiet = int(quiet)
        rc = int(rc)

        sel = selection
        ch = str(chain).strip()
        if ch and ch.lower() != "all":
            sel = "(%s) and chain %s" % (selection, ch)

        try:
            fasta = cmd.get_fastastr(sel)
        except Exception as e:
            raise CmdException("get_fastastr failed: %s" % e)

        blocks = DssrFunctions._parse_fastastr(fasta)
        if not blocks:
            raise CmdException('No FASTA sequence extracted from selection="%s"' % sel)

        out_lines = []
        for hdr, seq in blocks:
            s = "".join([c for c in seq.upper() if c.isalpha()])
            if rc:
                s = DssrFunctions._revcomp(s)

            if fmt == "fasta":
                out_lines.append(">" + hdr)
                out_lines.append(DssrFunctions._wrap_seq(s, wrap))
            else:
                out_lines.append(
                    hdr
                    + ": "
                    + (DssrFunctions._wrap_seq(s, wrap) if int(wrap) > 0 else s)
                )

        out = "\n".join(out_lines)
        if not quiet:
            print(out)
        return out

    @staticmethod
    def _restore_pretty_colors(objs):
        for o in objs:
            try:
                cmd.spectrum("resi", "rainbow", "(%s) and polymer" % o)
            except Exception:
                pass
            try:
                cmd.color("atomic", "(%s) and not polymer" % o)
            except Exception:
                pass

    @staticmethod
    def _clear_keep_molecules(apply_gray):
        """Clear only objects and selections created by this plugin.

        Older versions deleted every atomless PyMOL object, measurement, and
        named selection.  That could remove unrelated CGO/map objects and user
        work.  Molecular recoloring behavior is preserved, but deletion is
        now limited to names explicitly tracked by DSSR-PyMOL.
        """
        try:
            all_objs = list(cmd.get_object_list())
        except Exception:
            all_objs = []

        keep = []
        for o in all_objs:
            try:
                n = int(cmd.count_atoms(o))
            except Exception:
                n = 0
            if n > 0:
                keep.append(o)

        for o in list(_DSSR_BLOCK_OBJECTS):
            try:
                cmd.delete(o)
            except Exception:
                pass
        _DSSR_BLOCK_OBJECTS.clear()

        for selection_name in list(_DSSR_SELECTION_OBJECTS):
            try:
                cmd.delete(selection_name)
            except Exception:
                pass
        _DSSR_SELECTION_OBJECTS.clear()

        try:
            cmd.select("sele", "none")
        except Exception:
            pass

        if apply_gray:
            for o in keep:
                try:
                    cmd.color("gray", o)
                except Exception:
                    pass
        else:
            DssrFunctions._restore_pretty_colors(keep)

        try:
            cmd.zoom("all")
        except Exception:
            pass
        try:
            cmd.reset()
        except Exception:
            pass


class DssrGuiDialog(QtWidgets.QDialog if QtWidgets else object):
    def __init__(self):
        super(DssrGuiDialog, self).__init__()
        self.setWindowTitle("DSSR GUI")
        self.resize(1180, 760)

        self._current_feature = "pairs"
        self._cache_key = None
        self._cache_data = None
        self._items_all = []
        self._items_filtered = []
        self._page = 0
        self._query_conditions = []
        self._json_override = None
        self._json_override_context = None
        self._json_view_data = None
        self._json_text_cache = None
        self._json_box_data_id = None

        root = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QGridLayout()
        root.addLayout(top)

        self.obj_combo = QtWidgets.QComboBox()
        self.obj_combo.currentIndexChanged.connect(self._on_object_changed)

        self.state_combo = QtWidgets.QComboBox()
        self.state_combo.currentIndexChanged.connect(self._on_dssr_context_changed)
        self.count_btn = QtWidgets.QPushButton("count")
        self.count_btn.clicked.connect(self._count_states_clicked)

        self.exe_edit = QtWidgets.QLineEdit("x3dna-dssr")
        self.exe_edit.textChanged.connect(self._on_dssr_context_changed)

        self.color_edit = QtWidgets.QLineEdit("auto")
        self.name_edit = QtWidgets.QLineEdit("")

        self.precolor_cb = QtWidgets.QCheckBox("gray precolor")
        self.precolor_cb.setChecked(True)

        self.display_cb = QtWidgets.QCheckBox("display sticks")
        self.display_cb.setChecked(False)

        self.zoom_cb = QtWidgets.QCheckBox("zoom")
        self.zoom_cb.setChecked(True)

        self.showinfo_cb = QtWidgets.QCheckBox("show_info (pseudoknot)")
        self.showinfo_cb.setChecked(False)

        self.radius_spin = QtWidgets.QDoubleSpinBox()
        self.radius_spin.setMinimum(0.01)
        self.radius_spin.setMaximum(5.0)
        self.radius_spin.setSingleStep(0.05)
        self.radius_spin.setValue(0.25)

        self.status_label = QtWidgets.QLabel("")
        self.status_label.setWordWrap(True)

        self.block_file_combo = QtWidgets.QComboBox()
        self.block_file_combo.setEditable(True)
        self.block_file_combo.addItems(BLOCK_FEATURES)

        self.block_depth_spin = QtWidgets.QDoubleSpinBox()
        self.block_depth_spin.setMinimum(0.01)
        self.block_depth_spin.setMaximum(5.0)
        self.block_depth_spin.setSingleStep(0.05)
        self.block_depth_spin.setValue(0.5)

        self.make_blocks_btn = QtWidgets.QPushButton("make blocks")
        self.make_blocks_btn.clicked.connect(self._make_blocks_clicked)

        self.seq_btn = QtWidgets.QPushButton("seq view")
        self.seq_btn.clicked.connect(self._seq_view_clicked)

        top.addWidget(QtWidgets.QLabel("object"), 0, 0)
        top.addWidget(self.obj_combo, 0, 1, 1, 2)

        top.addWidget(QtWidgets.QLabel("state"), 0, 3)
        top.addWidget(self.state_combo, 0, 4)
        top.addWidget(self.count_btn, 0, 5)

        top.addWidget(QtWidgets.QLabel("exe"), 1, 0)
        top.addWidget(self.exe_edit, 1, 1, 1, 5)

        top.addWidget(QtWidgets.QLabel("color"), 2, 0)
        top.addWidget(self.color_edit, 2, 1)
        top.addWidget(QtWidgets.QLabel("name"), 2, 3)
        top.addWidget(self.name_edit, 2, 4)

        top.addWidget(QtWidgets.QLabel("block_file"), 3, 0)
        top.addWidget(self.block_file_combo, 3, 1)
        top.addWidget(QtWidgets.QLabel("block_depth"), 3, 3)
        top.addWidget(self.block_depth_spin, 3, 4)
        top.addWidget(self.make_blocks_btn, 3, 5, 1, 1)

        opts = QtWidgets.QHBoxLayout()
        opts.addWidget(self.precolor_cb)
        opts.addWidget(self.display_cb)
        opts.addWidget(self.zoom_cb)
        opts.addWidget(self.showinfo_cb)
        opts.addWidget(QtWidgets.QLabel("stick_radius"))
        opts.addWidget(self.radius_spin)
        opts.addStretch(1)
        root.addLayout(opts)

        root.addWidget(self.status_label)

        btn_area = QtWidgets.QWidget()
        btn_grid = QtWidgets.QGridLayout(btn_area)
        btn_grid.setContentsMargins(0, 0, 0, 0)
        btn_grid.setHorizontalSpacing(6)
        btn_grid.setVerticalSpacing(6)

        self._feature_buttons = {}
        cols = 6
        for idx, feat in enumerate(FEATURE_ORDER):
            b = QtWidgets.QPushButton(FEATURE_LABELS.get(feat, feat))
            b.setCheckable(True)
            b.clicked.connect(self._make_feature_handler(feat))
            r = idx // cols
            c = idx % cols
            btn_grid.addWidget(b, r, c)
            self._feature_buttons[feat] = b
        if "pairs" in self._feature_buttons:
            self._feature_buttons["pairs"].setChecked(True)

        root.addWidget(btn_area)

        bar = QtWidgets.QHBoxLayout()

        self.filter_edit = QtWidgets.QLineEdit("")
        self.filter_edit.setPlaceholderText("filter...")
        self.filter_edit.setMinimumWidth(170)
        self.filter_edit.textChanged.connect(self._on_filter_changed)

        self.page_size_spin = QtWidgets.QSpinBox()
        self.page_size_spin.setMinimum(50)
        self.page_size_spin.setMaximum(5000)
        self.page_size_spin.setValue(500)
        self.page_size_spin.valueChanged.connect(self._on_page_size_changed)

        self.prev_btn = QtWidgets.QPushButton("Prev")
        self.next_btn = QtWidgets.QPushButton("Next")
        self.prev_btn.clicked.connect(self._prev_page)
        self.next_btn.clicked.connect(self._next_page)

        self.page_label = QtWidgets.QLabel("")

        self.refresh_obj_btn = QtWidgets.QPushButton("refresh objects")
        self.refresh_obj_btn.clicked.connect(self.refresh_objects)
        self.refresh_list_btn = QtWidgets.QPushButton("refresh list")
        self.refresh_list_btn.clicked.connect(self._refresh_list_clicked)

        self.zoom_all_btn = QtWidgets.QPushButton("zoom all")
        self.zoom_all_btn.clicked.connect(self._zoom_all)

        self.reset_view_btn = QtWidgets.QPushButton("reset view")
        self.reset_view_btn.clicked.connect(self._reset_view)

        bar.addWidget(QtWidgets.QLabel("filter"))
        bar.addWidget(self.filter_edit, 2)
        bar.addWidget(QtWidgets.QLabel("max/page"))
        bar.addWidget(self.page_size_spin, 0)
        bar.addWidget(self.prev_btn)
        bar.addWidget(self.next_btn)
        bar.addWidget(self.page_label, 1)
        bar.addWidget(self.refresh_obj_btn)
        bar.addWidget(self.refresh_list_btn)
        bar.addWidget(self.seq_btn)
        bar.addWidget(self.zoom_all_btn)
        bar.addWidget(self.reset_view_btn)
        root.addLayout(bar)

        self.query_group = QtWidgets.QGroupBox("DSSR property query")
        query_outer = QtWidgets.QVBoxLayout(self.query_group)
        query_row = QtWidgets.QHBoxLayout()

        self.query_join_combo = QtWidgets.QComboBox()
        self.query_join_combo.addItems(["AND", "OR"])
        self.query_field_combo = QtWidgets.QComboBox()
        self.query_field_combo.setEditable(True)
        self.query_field_combo.setMinimumWidth(150)
        self.query_op_combo = QtWidgets.QComboBox()
        self.query_op_combo.addItems(
            ["=", "!=", "contains", "not contains", ">", ">=", "<", "<=", "exists", "not exists"]
        )
        self.query_value_edit = QtWidgets.QLineEdit()
        self.query_value_edit.setPlaceholderText("value, e.g. cWH, WC, 3.2, true")
        self.query_value_edit.returnPressed.connect(self._add_query_condition)

        self.query_add_btn = QtWidgets.QPushButton("add")
        self.query_add_btn.clicked.connect(self._add_query_condition)
        self.query_remove_btn = QtWidgets.QPushButton("undo condition")
        self.query_remove_btn.clicked.connect(self._remove_query_condition)
        self.query_clear_btn = QtWidgets.QPushButton("clear query")
        self.query_clear_btn.clicked.connect(self._clear_query_conditions)
        self.query_select_btn = QtWidgets.QPushButton("select matches in 3D")
        self.query_select_btn.clicked.connect(self._select_query_matches)

        query_row.addWidget(self.query_join_combo)
        query_row.addWidget(self.query_field_combo, 1)
        query_row.addWidget(self.query_op_combo)
        query_row.addWidget(self.query_value_edit, 2)
        query_row.addWidget(self.query_add_btn)
        query_row.addWidget(self.query_remove_btn)
        query_row.addWidget(self.query_clear_btn)
        query_row.addWidget(self.query_select_btn)
        query_outer.addLayout(query_row)
        self.query_summary_label = QtWidgets.QLabel("No property query — text filter remains available above")
        self.query_summary_label.setWordWrap(True)
        query_outer.addWidget(self.query_summary_label)
        root.addWidget(self.query_group)

        report_bar = QtWidgets.QHBoxLayout()
        self.report_btn = QtWidgets.QPushButton("Generate RNA Report")
        self.report_btn.clicked.connect(self._generate_rna_report_clicked)

        self.make_all_sel_btn = QtWidgets.QPushButton("make selections (all)")
        self.make_all_sel_btn.clicked.connect(self._make_all_selections_clicked)

        self.auto_color_btn = QtWidgets.QPushButton("auto color")
        self.auto_color_btn.clicked.connect(self._auto_color_clicked)

        self.clear_report_btn = QtWidgets.QPushButton("clear report")
        self.clear_report_btn.clicked.connect(self._clear_report)

        self.copy_json_btn = QtWidgets.QPushButton("copy JSON")
        self.copy_json_btn.clicked.connect(self._copy_json_clicked)
        self.save_json_btn = QtWidgets.QPushButton("save JSON")
        self.save_json_btn.clicked.connect(self._save_json_clicked)
        self.load_json_btn = QtWidgets.QPushButton("load JSON")
        self.load_json_btn.clicked.connect(self._load_json_clicked)

        report_bar.addWidget(self.report_btn)
        report_bar.addWidget(self.make_all_sel_btn)
        report_bar.addWidget(self.auto_color_btn)
        report_bar.addWidget(self.clear_report_btn)
        report_bar.addWidget(self.copy_json_btn)
        report_bar.addWidget(self.save_json_btn)
        report_bar.addWidget(self.load_json_btn)
        report_bar.addStretch(1)
        root.addLayout(report_bar)
        self.list_widget = QtWidgets.QListWidget()
        try:
            self.list_widget.setSelectionMode(
                QtWidgets.QAbstractItemView.ExtendedSelection
            )
        except Exception:
            pass
        self.list_widget.itemClicked.connect(self._on_item_clicked_preview)
        self.list_widget.itemDoubleClicked.connect(self._on_item_double_clicked)
        self.report_box = QtWidgets.QPlainTextEdit()
        self.report_box.setReadOnly(True)
        try:
            self.report_box.setPlaceholderText(
                "RNA Structure Summary will appear here..."
            )
        except Exception:
            pass

        self.details_box = QtWidgets.QPlainTextEdit()
        self.details_box.setReadOnly(True)
        self.details_box.setPlaceholderText(
            "Click a DSSR item to inspect every JSON property."
        )
        self.json_box = QtWidgets.QPlainTextEdit()
        self.json_box.setReadOnly(True)
        self.json_box.setPlaceholderText("Complete DSSR JSON will appear here.")
        self.data_tabs = QtWidgets.QTabWidget()
        self.data_tabs.addTab(self.report_box, "Report")
        self.data_tabs.addTab(self.details_box, "Item details")
        self.data_tabs.addTab(self.json_box, "Full JSON")
        self.data_tabs.currentChanged.connect(self._on_data_tab_changed)

        split = QtWidgets.QSplitter()
        split.setOrientation(QtCore.Qt.Horizontal)
        split.addWidget(self.list_widget)
        split.addWidget(self.data_tabs)
        try:
            split.setStretchFactor(0, 3)
            split.setStretchFactor(1, 2)
        except Exception:
            pass

        root.addWidget(split, 1)

        self.setStyleSheet(
            """
            QDialog { background: #ffffff; color: #263746; }
            QLabel, QCheckBox, QGroupBox { color: #263746; }
            QGroupBox {
                border: 1px solid #d4dde5; border-radius: 9px;
                margin-top: 8px; padding-top: 8px; background: #ffffff;
                font-weight: 600;
            }
            QGroupBox::title {
                subcontrol-origin: margin; left: 10px; padding: 0 6px;
                color: #24657d;
            }
            QPushButton, QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {
                color: #263746; background: #ffffff;
                border: 1px solid #bdcbd6; border-radius: 6px;
                padding: 4px 7px; min-height: 21px;
            }
            QPushButton:hover, QComboBox:hover, QLineEdit:focus {
                background: #eef9fd; border-color: #52b9da;
            }
            QPushButton:checked {
                color: #16475b; background: #d8f1fa; border-color: #55b9d9;
                font-weight: 700;
            }
            QPushButton:disabled, QComboBox:disabled {
                color: #98a6b1; background: #f3f5f7; border-color: #dce2e7;
            }
            QListWidget, QPlainTextEdit {
                color: #263746; background: #ffffff;
                border: 1px solid #d4dde5; border-radius: 7px;
                selection-color: #123d50; selection-background-color: #d6f0fa;
            }
            QComboBox QAbstractItemView {
                color: #263746; background: #ffffff;
                selection-color: #123d50; selection-background-color: #d6f0fa;
            }
            QTabWidget::pane { border: 1px solid #d4dde5; border-radius: 7px; }
            QTabBar::tab {
                color: #647581; background: #f2f5f7;
                border: 1px solid #d4dde5; padding: 6px 12px;
            }
            QTabBar::tab:selected { color: #174c61; background: #ffffff; }
            QToolTip { color: #263746; background: #ffffff;
                border: 1px solid #52b9da; }
            """
        )

        self.refresh_objects()
        self._update_state_combo()
        self.refresh_list()

    def _update_feature_buttons_state(self, dssr_data):
        """
        Dynamically enables/disables structural feature buttons based on
        the presence of features inside the loaded structure's DSSR output,
        and labels each button with its detected item count.
        """
        for feat, button in self._feature_buttons.items():
            label = FEATURE_LABELS.get(feat, feat)
            if not dssr_data:
                button.setEnabled(False)
                button.setText(label)
                continue

            if feat == "pseudoknot":
                count = ParsingAlgos._count_pseudoknot_layers(dssr_data)
            else:
                count = len(ParsingAlgos.feature_entries(dssr_data, feat))

            button.setEnabled(count > 0)
            button.setText("%s (%d)" % (label, count))

    def _enable_seq_view(self):
        try:
            cmd.set("seq_view", 1)
        except Exception:
            try:
                cmd.do("set seq_view, 1")
            except Exception:
                pass

    def _seq_view_clicked(self):
        self._enable_seq_view()
        sel = self._get_object_text()
        try:
            cmd.select("sele", "(%s)" % sel)
        except Exception:
            pass

    def _make_feature_handler(self, feat):
        def handler():
            self._current_feature = feat
            self._clear_query_conditions(render=False)
            for k, b in self._feature_buttons.items():
                b.setChecked(k == feat)
            self.refresh_list()

        return handler

    def _invalidate_dssr_cache(self):
        self._cache_key = None
        self._cache_data = None
        self._json_override = None
        self._json_override_context = None

    def _on_dssr_context_changed(self, *_args):
        self._invalidate_dssr_cache()

    def _on_object_changed(self, *_args):
        self._invalidate_dssr_cache()
        self._update_state_combo()

    def _refresh_list_clicked(self, *_args):
        self._invalidate_dssr_cache()
        self.refresh_list()

    def _count_states_clicked(self):
        obj = self._get_object_text()
        try:
            n = int(cmd.count_states(obj))
        except Exception:
            n = 0
        self._update_state_combo(force_n=n)
        try:
            QtWidgets.QMessageBox.information(
                self, "count_states", "count_states %s = %d" % (obj, n)
            )
        except Exception:
            pass

    def _get_object_text(self):
        txt = self.obj_combo.currentText().strip() if self.obj_combo.count() else ""
        return txt if txt else "all"

    def _update_state_combo(self, force_n=None):
        obj = self._get_object_text()
        try:
            n = int(force_n) if force_n is not None else int(cmd.count_states(obj))
        except Exception:
            n = 0

        current = self.state_combo.currentData() if self.state_combo.count() else None
        self.state_combo.clear()

        self.state_combo.addItem("current", -1)
        if n and n > 1:
            for s in range(1, n + 1):
                self.state_combo.addItem(str(s), s)

        if current is not None:
            idx = self.state_combo.findData(current)
            if idx >= 0:
                self.state_combo.setCurrentIndex(idx)

    def _get_state_value(self):
        data = self.state_combo.currentData()
        try:
            data = int(data)
        except Exception:
            data = -1
        if data == -1:
            try:
                return int(cmd.get_state())
            except Exception:
                return 1
        return data

    def refresh_objects(self):
        self._invalidate_dssr_cache()
        try:
            objs = cmd.get_object_list("enabled")
            if not objs:
                objs = cmd.get_object_list()
        except Exception:
            objs = []

        current = self.obj_combo.currentText().strip() if self.obj_combo.count() else ""
        self.obj_combo.clear()

        if not objs:
            self.obj_combo.addItem("all")
        else:
            for o in objs:
                self.obj_combo.addItem(o)

        if current:
            i = self.obj_combo.findText(current)
            if i >= 0:
                self.obj_combo.setCurrentIndex(i)
        else:
            try:
                enabled = cmd.get_object_list("enabled")
            except Exception:
                enabled = []
            if len(enabled) == 1:
                i = self.obj_combo.findText(enabled[0])
                if i >= 0:
                    self.obj_combo.setCurrentIndex(i)

        self._update_state_combo()

    def _zoom_all(self):
        try:
            cmd.zoom("all")
        except Exception:
            pass

    def _reset_view(self):
        apply_gray = True if self.precolor_cb.isChecked() else False
        DssrFunctions._clear_keep_molecules(apply_gray)
        self._invalidate_dssr_cache()
        self.refresh_objects()
        self.refresh_list()

    def _clear_report(self):
        try:
            self.report_box.setPlainText("")
        except Exception:
            pass

    @staticmethod
    def _json_pretty(value):
        return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)

    def _json_text_for_data(self, dssr_data):
        if dssr_data is self._json_view_data and self._json_text_cache is not None:
            return self._json_text_cache
        text = self._json_pretty(dssr_data)
        if dssr_data is self._json_view_data:
            self._json_text_cache = text
        return text

    def _render_json_view(self):
        data = self._json_view_data
        if data is None:
            return
        data_id = id(data)
        if self._json_box_data_id == data_id:
            return
        try:
            self.json_box.setPlainText(self._json_text_for_data(data))
            self._json_box_data_id = data_id
        except Exception as error:
            self.json_box.setPlainText("JSON rendering error: %s" % str(error))
            self._json_box_data_id = None

    def _on_data_tab_changed(self, _index):
        try:
            if self.data_tabs.currentWidget() is self.json_box:
                self._render_json_view()
        except Exception:
            pass

    def _refresh_json_view(self, dssr_data):
        if dssr_data is not self._json_view_data:
            self._json_view_data = dssr_data
            self._json_text_cache = None
            self._json_box_data_id = None
            try:
                self.json_box.setPlainText(
                    "Complete DSSR JSON is ready. Open this tab to render it."
                )
            except Exception:
                pass
        try:
            if self.data_tabs.currentWidget() is self.json_box:
                self._render_json_view()
        except Exception:
            pass

    def _copy_json_clicked(self):
        try:
            sel_obj, exe, st, precolor_on = self._get_dssr_context()
            data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            text = self._json_text_for_data(data)
            QtWidgets.QApplication.clipboard().setText(text)
            self.status_label.setText("Complete DSSR JSON copied to clipboard")
        except Exception as error:
            self.status_label.setText("copy JSON error: %s" % str(error))

    def _save_json_clicked(self):
        try:
            path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save DSSR JSON",
                "dssr-analysis.json",
                "JSON files (*.json);;All files (*)",
            )
            if not path:
                return
            sel_obj, exe, st, precolor_on = self._get_dssr_context()
            data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(self._json_text_for_data(data))
                handle.write("\n")
            self.status_label.setText("DSSR JSON saved: %s" % path)
        except Exception as error:
            QtWidgets.QMessageBox.critical(self, "Save DSSR JSON", str(error))

    def _load_json_clicked(self):
        try:
            path, _selected_filter = QtWidgets.QFileDialog.getOpenFileName(
                self,
                "Load DSSR JSON",
                "",
                "JSON files (*.json);;All files (*)",
            )
            if not path:
                return
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            if not isinstance(data, dict):
                raise CmdException("DSSR JSON root must be an object")
            context = (self._get_object_text(), self._get_state_value())
            self._cache_key = None
            self._cache_data = data
            self._json_override = data
            self._json_override_context = context
            self._refresh_json_view(data)
            self.refresh_list()
            self.status_label.setText(
                "Loaded DSSR JSON for %s state %s: %s"
                % (context[0], str(context[1]), path)
            )
        except Exception as error:
            QtWidgets.QMessageBox.critical(self, "Load DSSR JSON", str(error))

    @staticmethod
    def _query_scalar_paths(value, prefix="", depth=0):
        paths = set()
        if not isinstance(value, dict):
            return paths
        for key, child in value.items():
            key = str(key)
            path = "%s.%s" % (prefix, key) if prefix else key
            paths.add(path)
            if isinstance(child, dict) and depth < 2:
                paths |= DssrGuiDialog._query_scalar_paths(child, path, depth + 1)
        return paths

    @staticmethod
    def _query_get_path(entry, path):
        marker = object()
        value = entry
        for part in str(path).split("."):
            if not isinstance(value, dict) or part not in value:
                return marker, False
            value = value[part]
        return value, True

    @staticmethod
    def _query_coerce(value):
        if not isinstance(value, str):
            return value
        text = value.strip()
        if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
            text = text[1:-1]
        lowered = text.lower()
        if lowered == "true":
            return True
        if lowered == "false":
            return False
        if lowered in ("none", "null"):
            return None
        try:
            return float(text)
        except Exception:
            return text

    @staticmethod
    def _query_compare(entry, field, operator, raw_expected):
        actual, exists = DssrGuiDialog._query_get_path(entry, field)
        operator = str(operator).lower().strip()
        if operator == "exists":
            return exists and actual is not None
        if operator == "not exists":
            return (not exists) or actual is None
        if not exists:
            return False

        expected = DssrGuiDialog._query_coerce(raw_expected)
        actual_coerced = DssrGuiDialog._query_coerce(actual)
        if isinstance(actual, (dict, list, tuple)):
            actual_text = json.dumps(actual, sort_keys=True, ensure_ascii=False)
        else:
            actual_text = str(actual)
        expected_text = str(expected)

        if operator == "contains":
            return expected_text.lower() in actual_text.lower()
        if operator == "not contains":
            return expected_text.lower() not in actual_text.lower()
        if operator in ("=", "!="):
            if isinstance(actual_coerced, (int, float, bool)) or actual_coerced is None:
                matched = actual_coerced == expected
            else:
                matched = str(actual_coerced).lower() == expected_text.lower()
            return matched if operator == "=" else not matched

        try:
            left = float(actual_coerced)
            right = float(expected)
        except Exception:
            left = str(actual_coerced).lower()
            right = expected_text.lower()
        if operator == ">":
            return left > right
        if operator == ">=":
            return left >= right
        if operator == "<":
            return left < right
        if operator == "<=":
            return left <= right
        return False

    def _entry_for_feature_index(self, dssr_data, feature, index):
        feature = str(feature).strip()
        index = int(index)
        if feature == "pseudoknot":
            dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
            layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket)
            keys = sorted(layers.keys())
            if 1 <= index <= len(keys):
                key = keys[index - 1]
                return {
                    "layer": key,
                    "pair_count": len(layers[key]),
                    "pairs": layers[key],
                }
            return None
        entries = ParsingAlgos.feature_entries(dssr_data, feature)
        if 1 <= index <= len(entries):
            entry = entries[index - 1]
            return entry if isinstance(entry, dict) else {"value": entry}
        return None

    def _update_query_fields(self, dssr_data):
        current = self.query_field_combo.currentText().strip()
        fields = set()
        for index, _text in self._items_all[:500]:
            entry = self._entry_for_feature_index(
                dssr_data, self._current_feature, index
            )
            fields |= self._query_scalar_paths(entry)
        self.query_field_combo.blockSignals(True)
        self.query_field_combo.clear()
        self.query_field_combo.addItems(sorted(fields, key=lambda item: item.lower()))
        if current and current in fields:
            self.query_field_combo.setCurrentText(current)
        self.query_field_combo.blockSignals(False)
        enabled = bool(fields)
        for widget in (
            self.query_field_combo,
            self.query_op_combo,
            self.query_value_edit,
            self.query_add_btn,
            self.query_select_btn,
        ):
            widget.setEnabled(enabled)

    def _query_summary_text(self):
        if not self._query_conditions:
            return "No property query — text filter remains available above"
        chunks = []
        for pos, condition in enumerate(self._query_conditions):
            join, field, operator, value = condition
            prefix = "" if pos == 0 else (join + " ")
            if operator in ("exists", "not exists"):
                chunks.append("%s%s %s" % (prefix, field, operator))
            else:
                chunks.append("%s%s %s %s" % (prefix, field, operator, value))
        return "Query: " + " ".join(chunks)

    def _add_query_condition(self):
        field = self.query_field_combo.currentText().strip()
        operator = self.query_op_combo.currentText().strip()
        value = self.query_value_edit.text().strip()
        if not field:
            return
        if operator not in ("exists", "not exists") and value == "":
            self.status_label.setText("Enter a query value")
            return
        join = self.query_join_combo.currentText().strip().upper() or "AND"
        self._query_conditions.append((join, field, operator, value))
        self.query_value_edit.clear()
        self.query_summary_label.setText(self._query_summary_text())
        self._page = 0
        self._render_list()

    def _remove_query_condition(self):
        if self._query_conditions:
            self._query_conditions.pop()
        self.query_summary_label.setText(self._query_summary_text())
        self._page = 0
        self._render_list()

    def _clear_query_conditions(self, render=True):
        self._query_conditions = []
        try:
            self.query_summary_label.setText(self._query_summary_text())
        except Exception:
            pass
        if render:
            self._page = 0
            self._render_list()

    def _entry_matches_query(self, entry):
        if not self._query_conditions:
            return True
        result = None
        for join, field, operator, value in self._query_conditions:
            matched = self._query_compare(entry, field, operator, value)
            if result is None:
                result = matched
            elif join == "OR":
                result = bool(result or matched)
            else:
                result = bool(result and matched)
        return bool(result)

    def _select_query_matches(self):
        try:
            self._render_list()
            if not self._items_filtered:
                raise CmdException("No items match the current text/property query")
            sel_obj, exe, st, precolor_on = self._get_dssr_context()
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            parts = []
            for index, _text in self._items_filtered:
                try:
                    parts.append(
                        "(%s)"
                        % ParsingAlgos._build_residue_sel_from_dssr(
                            dssr_data,
                            self._current_feature,
                            index,
                        )
                    )
                except Exception:
                    pass
            if not parts:
                raise CmdException("Matched JSON items could not be mapped to PyMOL")
            name = "dssr_query_matches"
            cmd.select(name, "((%s) and (%s))" % (sel_obj, " or ".join(parts)))
            _DSSR_SELECTION_OBJECTS.add(name)
            cmd.color("cyan", name)
            if self.zoom_cb.isChecked():
                cmd.zoom(name, buffer=4.0)
            self.status_label.setText(
                "Selected %d matching DSSR item(s) as %s"
                % (len(self._items_filtered), name)
            )
        except Exception as error:
            self.status_label.setText("query selection error: %s" % str(error))

    def _show_item_details(self, dssr_data, feature, index):
        entry = self._entry_for_feature_index(dssr_data, feature, index)
        payload = {
            "feature": feature,
            "array_index": int(index),
            "data": entry,
        }
        self.details_box.setPlainText(self._json_pretty(payload))
        self.data_tabs.setCurrentWidget(self.details_box)

    def _get_dssr_context(self):
        sel_obj = self._get_object_text()
        exe = self.exe_edit.text().strip() or "x3dna-dssr"
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0
        return sel_obj, exe, st, precolor_on

    def _append_report(self, text):
        try:
            cur = self.report_box.toPlainText()
        except Exception:
            cur = ""
        if cur:
            out = cur.rstrip("\n") + "\n\n" + str(text).rstrip("\n") + "\n"
        else:
            out = str(text).rstrip("\n") + "\n"
        try:
            self.report_box.setPlainText(out)
        except Exception:
            pass

    def _generate_rna_report_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            report = ParsingAlgos._format_rna_summary_text(dssr_data)
            try:
                self.report_box.setPlainText(report + "\n")
            except Exception:
                self._append_report(report)
        except Exception as e:
            msg = "RNA report error: %s" % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, "RNA Report", msg)
            except Exception:
                pass

    def _make_all_selections_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            made, skipped = self._make_all_selections(dssr_data, sel_obj)
            self._append_report(
                "Created selections: %s" % (" ".join(made) if made else "(none)")
            )
            if skipped:
                self._append_report("Skipped (empty): %s" % (" ".join(skipped)))
        except Exception as e:
            msg = "make selections error: %s" % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, "make selections", msg)
            except Exception:
                pass

    def _auto_color_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            self._make_all_selections(dssr_data, sel_obj)
            self._apply_auto_colors()
            self._append_report(
                "Auto color applied: stems green, hairpins blue, pseudoknots red, aminors purple"
            )
        except Exception as e:
            msg = "auto color error: %s" % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, "auto color", msg)
            except Exception:
                pass

    def _make_all_selections(self, dssr_data, obj_sel):
        targets = [
            ("pairs", "pairs_all"),
            ("hairpins", "hairpins_all"),
            ("stems", "stems_all"),
            ("bulges", "bulges_all"),
            ("junctions", "junctions_all"),
            ("pseudoknot", "pseudoknots_all"),
            ("aminors", "aminors_all"),
            ("stacks", "stacks_all"),
            ("uturns", "uturns_all"),
        ]

        made = []
        skipped = []

        for feat, name in targets:
            residues = ParsingAlgos._collect_residues_all(dssr_data, feat)
            sel_core = ParsingAlgos._compact_sel_from_residues(residues)
            if not sel_core:
                skipped.append(name)
                try:
                    cmd.delete(name)
                except Exception:
                    pass
                continue

            expr = "((%s) and (%s))" % (obj_sel, sel_core)
            try:
                cmd.select(name, expr)
                _DSSR_SELECTION_OBJECTS.add(str(name))
                made.append(name)
            except Exception:
                skipped.append(name)
                try:
                    cmd.delete(name)
                except Exception:
                    pass

        return made, skipped

    def _apply_auto_colors(self):
        try:
            cmd.color("green", "stems_all")
        except Exception:
            pass
        try:
            cmd.color("blue", "hairpins_all")
        except Exception:
            pass
        try:
            cmd.color("red", "pseudoknots_all")
        except Exception:
            pass
        try:
            cmd.color("purple", "aminors_all")
        except Exception:
            pass

    def _big_object_warning(self, sel):
        thresh = 50000
        try:
            n_atoms = int(cmd.count_atoms(sel))
        except Exception:
            n_atoms = 0
        if n_atoms >= thresh:
            self.status_label.setText(
                "Warning: Large selection (%d atoms). DSSR analysis may be slow. Consider selecting a specific chain."
                % n_atoms
            )
        else:
            self.status_label.setText("")

    def _get_dssr_data(self, selection, state, exe, precolor_on):
        override_context = (str(selection), int(state))
        if (
            self._json_override is not None
            and self._json_override_context == override_context
        ):
            if int(precolor_on):
                try:
                    cmd.color("gray", selection)
                except Exception:
                    pass
            self._refresh_json_view(self._json_override)
            return self._json_override

        try:
            atom_count = int(cmd.count_atoms(selection, state=state))
        except Exception:
            try:
                atom_count = int(cmd.count_atoms(selection))
            except Exception:
                atom_count = -1
        try:
            extent = cmd.get_extent(selection, state=state)
            extent_key = tuple(
                round(float(value), 3)
                for point in extent
                for value in point
            )
        except Exception:
            extent_key = ()

        cache_key = (
            str(selection),
            int(state),
            str(exe),
            atom_count,
            extent_key,
        )
        if self._cache_key == cache_key and self._cache_data is not None:
            if int(precolor_on):
                try:
                    cmd.color("gray", selection)
                except Exception:
                    pass
            self._refresh_json_view(self._cache_data)
            return self._cache_data

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmpfilepdb = tmp.name
        tmp.close()

        try:
            cmd.save(tmpfilepdb, selection, state)
            if int(precolor_on):
                cmd.color("gray", selection)
            data = HelperFunctions.run_dssr_json(tmpfilepdb, exe)
        finally:
            try:
                os.remove(tmpfilepdb)
            except OSError:
                pass

        self._cache_key = cache_key
        self._cache_data = data
        self._refresh_json_view(data)
        return data

    def _on_filter_changed(self, _):
        self._page = 0
        self._render_list()

    def _on_page_size_changed(self, _):
        self._page = 0
        self._render_list()

    def _prev_page(self):
        if self._page > 0:
            self._page -= 1
            self._render_list()

    def _next_page(self):
        pages = self._total_pages()
        if self._page + 1 < pages:
            self._page += 1
            self._render_list()

    def _total_pages(self):
        page_size = int(self.page_size_spin.value())
        total = len(self._items_filtered)
        if page_size <= 0:
            return 1
        pages = (total + page_size - 1) // page_size
        return max(1, pages)

    def refresh_list(self):
        sel = self._get_object_text()
        feat = self._current_feature
        exe = self.exe_edit.text().strip() or "x3dna-dssr"
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0

        self._big_object_warning(sel)

        self.list_widget.clear()

        # Guard Check: If there is no structure loaded in PyMOL, do not trigger DSSR or update with loading status.
        if not cmd.get_object_list():
            msg = "No structure loaded. Please load a PDB/CIF file before running DSSR-PyMOL"
            self.list_widget.addItem(
                "Please load a PDB/CIF file before running DSSR-PyMOL"
            )
            self.status_label.setText(msg)
            self._update_feature_buttons_state(None)
            return

        self.list_widget.addItem("loading...")
        QtWidgets.QApplication.processEvents()

        try:
            dssr_data = self._get_dssr_data(sel, st, exe, precolor_on)
            self._update_feature_buttons_state(dssr_data)

            items = []

            if feat == "pseudoknot":
                dotbracket = ParsingAlgos._extract_dotbracket(dssr_data)
                nts_list = dssr_data.get("nts", None)
                if nts_list is None:
                    self.list_widget.clear()
                    self.list_widget.addItem("No nucleotides found in DSSR output.")
                    return

                layers = ParsingAlgos.parse_dotbracket_pseudoknots(dotbracket)
                if not layers:
                    self.list_widget.clear()
                    self.list_widget.addItem(
                        "No pseudoknot layers found in this structure."
                    )
                    return

                layer_keys = sorted(layers.keys())
                for j, k in enumerate(layer_keys, 1):
                    line = "%d: layer key=%s pairs=%d" % (j, str(k), len(layers[k]))
                    items.append((j, line))

            else:
                if feat not in FEATURE_MAP:
                    raise CmdException('Unknown feature "%s"' % feat)
                json_key = FEATURE_MAP[feat]
                feature_list = ParsingAlgos.feature_entries(dssr_data, feat)
                if not feature_list:
                    nice_names = {
                        "pairs": "base pairs",
                        "stems": "stems",
                        "helices": "helices",
                        "stacks": "stacks",
                        "nonstack": "non-stacking nucleotides",
                        "coaxstacks": "coaxial stacks",
                        "atom2bases": "atom-to-base interactions",
                        "aminors": "A-minor interactions",
                        "splayunits": "splayed units",
                        "hairpins": "hairpin loops",
                        "bulges": "bulge loops",
                        "iloops": "internal loops",
                        "internal": "internal loops",
                        "junctions": "junction loops",
                        "sssegments": "single-stranded segments",
                        "ssSegments": "single-stranded segments",
                        "multiplets": "multiplets",
                        "nts": "nucleotides",
                        "pseudoknot": "pseudoknot layers",
                        "uturns": "U-turns",
                    }
                    feat_name = nice_names.get(feat, feat)
                    self.list_widget.clear()
                    self.list_widget.addItem(
                        "No %s found in this structure." % feat_name
                    )
                    self.page_label.setText("items 0, filtered 0, page 1/1")
                    self.prev_btn.setEnabled(False)
                    self.next_btn.setEnabled(False)
                    return

                total = len(feature_list)
                for i in range(total):
                    line = ParsingAlgos._preview_entry(feat, feature_list[i], i + 1)
                    items.append((i + 1, line))

            self._items_all = items
            self._update_query_fields(dssr_data)
            self._page = 0
            self._render_list()

        except Exception as e:
            self._items_all = []
            self._items_filtered = []
            self.list_widget.clear()
            self._update_feature_buttons_state(None)
            self.list_widget.addItem("ERROR: %s" % str(e))
            try:
                print("dssr_gui error: %s" % str(e))
            except Exception:
                pass

    def _render_list(self):
        q = self.filter_edit.text().strip().lower()
        filtered = list(self._items_all)
        if self._query_conditions:
            dssr_data = self._json_override or self._cache_data or {}
            query_filtered = []
            for item in filtered:
                entry = self._entry_for_feature_index(
                    dssr_data, self._current_feature, item[0]
                )
                if entry is not None and self._entry_matches_query(entry):
                    query_filtered.append(item)
            filtered = query_filtered
        if q:
            filtered = [it for it in filtered if q in it[1].lower()]
        self._items_filtered = filtered

        total_all = len(self._items_all)
        total_f = len(self._items_filtered)

        page_size = int(self.page_size_spin.value())
        pages = self._total_pages()

        if self._page >= pages:
            self._page = max(0, pages - 1)

        start = self._page * page_size
        end = start + page_size
        page_items = self._items_filtered[start:end]

        self.list_widget.clear()
        for idx, text in page_items:
            it = QtWidgets.QListWidgetItem(text)
            it.setData(QtCore.Qt.UserRole, int(idx))
            self.list_widget.addItem(it)

        self.page_label.setText(
            "items %d, filtered %d, page %d/%d"
            % (total_all, total_f, self._page + 1, pages)
        )
        self.prev_btn.setEnabled(self._page > 0)
        self.next_btn.setEnabled(self._page + 1 < pages)

        try:
            self.setWindowTitle("DSSR GUI - %s" % self._current_feature)
        except Exception:
            pass

    def _on_item_clicked_preview(self, item):
        self._enable_seq_view()

        data = item.data(QtCore.Qt.UserRole)
        if data is None:
            return
        try:
            idx = int(data)
        except Exception:
            return

        sel_obj = self._get_object_text()
        feat = self._current_feature
        exe = self.exe_edit.text().strip() or "x3dna-dssr"
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0

        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            sel_str = ParsingAlgos._build_residue_sel_from_dssr(
                dssr_data, feat, idx
            )
            cmd.select("sele", "((%s) and (%s))" % (sel_obj, sel_str))
            cmd.select("sele", "byres (sele)")
            self._show_item_details(dssr_data, feat, idx)
        except Exception as e:
            try:
                self.status_label.setText("preview error: %s" % str(e))
            except Exception:
                pass

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

        nm = self.name_edit.text().strip()
        if not nm:
            nm = "%s%d" % (feat, idx)

        col = self.color_edit.text().strip() or "auto"
        precolor_on = 1 if self.precolor_cb.isChecked() else 0
        display_on = 1 if self.display_cb.isChecked() else 0
        zoom_on = 1 if self.zoom_cb.isChecked() else 0
        showinfo_on = 1 if self.showinfo_cb.isChecked() else 0
        radius = float(self.radius_spin.value())

        try:
            DssrFunctions.dssr(
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
        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, "DSSR GUI error", str(e))
            except Exception:
                pass
            try:
                print("dssr_gui select error: %s" % str(e))
            except Exception:
                pass

    def _make_blocks_clicked(self):
        self._enable_seq_view()

        sel_obj = self._get_object_text()
        exe = self.exe_edit.text().strip() or "x3dna-dssr"
        st = self._get_state_value()

        block_file = self.block_file_combo.currentText().strip() or "face"
        block_depth = float(self.block_depth_spin.value())

        base = self.name_edit.text().strip()
        if not base:
            base = "blk"

        items = []
        try:
            items = list(self.list_widget.selectedItems())
        except Exception:
            items = []

        dssr_data = None
        if items:
            try:
                precolor_on = 1 if self.precolor_cb.isChecked() else 0
                dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            except Exception as e:
                try:
                    QtWidgets.QMessageBox.critical(self, "make blocks error", str(e))
                except Exception:
                    pass
                return

        made_names = []

        try:
            if items:
                feat = self._current_feature
                for it in items:
                    data = it.data(QtCore.Qt.UserRole)
                    if data is None:
                        continue
                    try:
                        idx = int(data)
                    except Exception:
                        continue

                    sel_str = ParsingAlgos._build_residue_sel_from_dssr(
                        dssr_data, feat, idx
                    )
                    sel_for_block = "byres (%s)" % sel_str

                    obj_name = "%s_%s_%d" % (base, feat, idx)
                    try:
                        cmd.delete(obj_name)
                    except Exception:
                        pass

                    DssrFunctions.dssr_block(
                        selection=sel_for_block,
                        state=st,
                        block_file=block_file,
                        block_depth=block_depth,
                        name=obj_name,
                        exe=exe,
                        quiet=1,
                    )

                    _DSSR_BLOCK_OBJECTS.add(obj_name)
                    made_names.append(obj_name)

            else:
                try:
                    n = int(cmd.count_atoms("sele"))
                except Exception:
                    n = 0
                if n > 0:
                    sel_for_block = "byres (sele)"
                else:
                    sel_for_block = "(%s)" % sel_obj

                obj_name = "%s_sel" % base
                try:
                    cmd.delete(obj_name)
                except Exception:
                    pass

                DssrFunctions.dssr_block(
                    selection=sel_for_block,
                    state=st,
                    block_file=block_file,
                    block_depth=block_depth,
                    name=obj_name,
                    exe=exe,
                    quiet=1,
                )

                _DSSR_BLOCK_OBJECTS.add(obj_name)
                made_names.append(obj_name)

            if made_names and self.zoom_cb.isChecked():
                try:
                    cmd.zoom("(%s)" % " or ".join(made_names))
                except Exception:
                    try:
                        cmd.zoom("all")
                    except Exception:
                        pass

        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, "make blocks error", str(e))
            except Exception:
                pass
            try:
                print("make blocks error: %s" % str(e))
            except Exception:
                pass

    @staticmethod
    def dssr_gui():
        global _DSSR_GUI_DIALOG
        if QtWidgets is None or QtCore is None:
            raise CmdException("Qt is not available in this PyMOL build")

        # Check if a structure is available before showing GUI or running DSSR
        no_struct = not cmd.get_object_list()

        if _DSSR_GUI_DIALOG is None:
            _DSSR_GUI_DIALOG = DssrGuiDialog()
        _DSSR_GUI_DIALOG.show()
        _DSSR_GUI_DIALOG.raise_()
        _DSSR_GUI_DIALOG.activateWindow()

        # Visual pop-up dialog warning and console diagnostic print if no structure is available
        if no_struct:
            msg = "No structure loaded. Please load a PDB/CIF file before running DSSR-PyMOL"
            print(msg)
            _DSSR_GUI_DIALOG.status_label.setText(msg)
            QtWidgets.QMessageBox.warning(_DSSR_GUI_DIALOG, "No Structure Loaded", msg)

dssr_select = DssrFunctions.dssr_select
dssr_gui = DssrGuiDialog.dssr_gui
dssr_block = DssrFunctions.dssr_block
dssr_seq = DssrFunctions.dssr_seq

cmd.extend("dssr_select", dssr_select)
cmd.extend("dssr_gui", dssr_gui)
cmd.extend("dssr_block", dssr_block)

# Restore tab-completion of arguments for dssr_block
try:
    cmd.auto_arg[0].update(
        {
            "dssr_block": cmd.auto_arg[0]["zoom"],
        }
    )
    cmd.auto_arg[2].update(
        {
            "dssr_block": [
                cmd.Shortcut(BLOCK_FEATURES),
                "block_file",
                "",
            ],
        }
    )
except Exception:
    pass



def __init_plugin__(app=None):
    addmenuitemqt("DSSR", dssr_gui)

try:
    print("Loaded DSSR helper %s from: %s" % (__DSSR_PLUGIN_VERSION__, __file__))
except Exception:
    pass

# ============================================================================
# Pure-Python RNA 2D extension for DSSR-PyMOL
#
# This extension is intentionally appended to the original plugin without
# removing or rewriting any of the original source above.  It replaces the
# external Java/Jmol/VARNA runtime requirement with an in-process Python/Qt
# viewer.  The implementation is original Python code.  It follows established
# RNA-visualization architecture (model -> layout -> scene -> interaction), but
# does not copy or translate VARNA's GPL source code.
# ============================================================================

import math
import random

try:
    from pymol.Qt import QtGui
except Exception:
    QtGui = None

try:
    from pymol.Qt import QtSvg
except Exception:
    QtSvg = None

__DSSR_PY2D_EXTENSION_VERSION__ = "v1.0.0-pure-python"
_DSSR_PY2D_DIALOGS = []


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
            # Prefer a conventional one-letter base when one is present.
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
            chain, resi = ParsingAlgos.parse_nt_id(nt_id)
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
                    ParsingAlgos._extract_dotbracket(dssr_data),
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
                seq_no_breaks = (seq_no_breaks + derived + ("N" * target_len))[:target_len]
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
            src = nts_json[i] if i < len(nts_json) and isinstance(nts_json[i], dict) else {}
            nt_id = str(src.get("nt_id", ""))
            chain, resi = cls._safe_parse_nt_id(nt_id)
            base = model.sequence[i] if i < len(model.sequence) else cls._base_from_nt_entry(src)
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
                        "Unmatched closing bracket %s at nucleotide %d" % (char, idx + 1)
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
                    "Unmatched opening symbol %s at nucleotide %d"
                    % (opener, idx + 1)
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
        a pair already in the scaffold.  All excluded pairs are still rendered
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
    def _graph_edges(model):
        edges = []
        n = len(model.nts)
        for i in range(n - 1):
            if i not in model.chain_breaks:
                edges.append((i, i + 1, 44.0, 1.00))
        for pair in model.planar_secondary_pairs():
            layer = int(pair.get("layer", 0))
            weight = 1.35 if layer == 0 else 0.95
            target = 39.0 if layer == 0 else 52.0
            edges.append((int(pair["i"]), int(pair["j"]), target, weight))
        return edges

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
    def radiate(model):
        """Draw a planar stem/loop scaffold and overlay pseudoknots later.

        This is an original, dependency-free radial layout.  Consecutive base
        pairs become ladder-like stems, child stems fan out from multiloops, and
        unpaired hairpin residues follow a long circular arc.  It deliberately
        keeps pseudoknot and tertiary edges out of the scaffold so they can be
        drawn on top without destroying the main secondary-structure topology.
        """
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
                if math.hypot(*average_direction) < 0.1:
                    chord = (b[0] - a[0], b[1] - a[1])
                    average_direction = (-chord[1], chord[0])
                average_direction = Dssr2DLayout._v_norm(average_direction)

                # A run directly enclosed by one base pair is a hairpin loop.
                # Use the long arc of a circle so bases remain well separated.
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
                        Dssr2DLayout._v_mul(
                            average_direction, center_distance
                        ),
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

                # Generic loop or linker: a smooth quadratic arc between anchors.
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

        # Process each chain independently so no artificial backbone connection
        # is introduced across an '&' separator.
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

        # Defensive fallback for unusual disconnected annotations.
        fallback = Dssr2DLayout.circular(model)
        for i in range(n):
            if positions[i] is None:
                positions[i] = fallback[i]

        center_x = sum(point[0] for point in positions) / float(n)
        center_y = sum(point[1] for point in positions) / float(n)
        return [
            (point[0] - center_x, point[1] - center_y) for point in positions
        ]

    @staticmethod
    def force(model):
        n = len(model.nts)
        if n <= 2:
            return Dssr2DLayout.linear(model)
        if n > 900:
            # A quadratic force solver is intentionally avoided for very large
            # structures.  Circular remains deterministic and immediately usable.
            return Dssr2DLayout.circular(model)

        initial = Dssr2DLayout.circular(model)
        rng = random.Random(1731 + n)
        pos = [
            [x + rng.uniform(-3.0, 3.0), y + rng.uniform(-3.0, 3.0)]
            for x, y in initial
        ]
        vel = [[0.0, 0.0] for _ in range(n)]
        edges = Dssr2DLayout._graph_edges(model)

        if n <= 120:
            iterations = 260
        elif n <= 300:
            iterations = 180
        elif n <= 600:
            iterations = 100
        else:
            iterations = 55

        repulsion = 18000.0
        max_repulsion_distance = 240.0
        max_repulsion_sq = max_repulsion_distance * max_repulsion_distance
        spring_k = 0.045
        gravity = 0.0018

        for iteration in range(iterations):
            forces = [[-gravity * p[0], -gravity * p[1]] for p in pos]

            # Pairwise repulsion.  The cutoff keeps medium-sized RNAs responsive.
            for i in range(n):
                xi, yi = pos[i]
                for j in range(i + 1, n):
                    dx = xi - pos[j][0]
                    dy = yi - pos[j][1]
                    d2 = dx * dx + dy * dy
                    if d2 > max_repulsion_sq:
                        continue
                    if d2 < 4.0:
                        dx += rng.uniform(-1.0, 1.0)
                        dy += rng.uniform(-1.0, 1.0)
                        d2 = max(4.0, dx * dx + dy * dy)
                    d = math.sqrt(d2)
                    f = repulsion / d2
                    fx = f * dx / d
                    fy = f * dy / d
                    forces[i][0] += fx
                    forces[i][1] += fy
                    forces[j][0] -= fx
                    forces[j][1] -= fy

            # Backbone and base-pair springs.
            for i, j, target, weight in edges:
                dx = pos[j][0] - pos[i][0]
                dy = pos[j][1] - pos[i][1]
                d2 = dx * dx + dy * dy
                if d2 < 1.0e-8:
                    continue
                d = math.sqrt(d2)
                f = spring_k * weight * (d - target)
                fx = f * dx / d
                fy = f * dy / d
                forces[i][0] += fx
                forces[i][1] += fy
                forces[j][0] -= fx
                forces[j][1] -= fy

            # Light local straightening keeps stems/backbones readable without
            # imposing a rigid template.
            for i in range(1, n - 1):
                if (i - 1) in model.chain_breaks or i in model.chain_breaks:
                    continue
                mx = (pos[i - 1][0] + pos[i + 1][0]) * 0.5
                my = (pos[i - 1][1] + pos[i + 1][1]) * 0.5
                forces[i][0] += (mx - pos[i][0]) * 0.012
                forces[i][1] += (my - pos[i][1]) * 0.012

            cooling = 1.0 - 0.70 * (float(iteration) / max(1.0, iterations - 1.0))
            max_step = 9.0 * cooling + 1.5
            for i in range(n):
                vel[i][0] = (vel[i][0] + forces[i][0]) * 0.78
                vel[i][1] = (vel[i][1] + forces[i][1]) * 0.78
                speed = math.hypot(vel[i][0], vel[i][1])
                if speed > max_step:
                    scale = max_step / speed
                    vel[i][0] *= scale
                    vel[i][1] *= scale
                pos[i][0] += vel[i][0]
                pos[i][1] += vel[i][1]

        # Center the final coordinates.
        cx = sum(p[0] for p in pos) / n
        cy = sum(p[1] for p in pos) / n
        return [(p[0] - cx, p[1] - cy) for p in pos]

    @staticmethod
    def compute(model, algorithm):
        name = str(algorithm or "radiate").strip().lower()
        if name in ("linear", "line", "feynman"):
            return Dssr2DLayout.linear(model)
        if name in ("circular", "circle"):
            return Dssr2DLayout.circular(model)
        if name in ("force", "spring", "graph"):
            return Dssr2DLayout.force(model)
        return Dssr2DLayout.radiate(model)


if QtWidgets is not None and QtCore is not None and QtGui is not None:

    class Dssr2DGraphicsView(QtWidgets.QGraphicsView):
        def __init__(self, scene, parent=None):
            super(Dssr2DGraphicsView, self).__init__(scene, parent)
            try:
                self.setRenderHints(
                    QtGui.QPainter.Antialiasing
                    | QtGui.QPainter.TextAntialiasing
                    | QtGui.QPainter.SmoothPixmapTransform
                )
            except Exception:
                try:
                    self.setRenderHint(QtGui.QPainter.Antialiasing, True)
                except Exception:
                    pass
            try:
                self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
                self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
                self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorViewCenter)
            except Exception:
                pass
            try:
                self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("white")))
            except Exception:
                pass

        def wheelEvent(self, event):
            try:
                delta = event.angleDelta().y()
            except Exception:
                try:
                    delta = event.delta()
                except Exception:
                    delta = 0
            factor = 1.18 if delta > 0 else (1.0 / 1.18)
            try:
                self.scale(factor, factor)
                event.accept()
            except Exception:
                super(Dssr2DGraphicsView, self).wheelEvent(event)


    class Dssr2DEdgeItem(QtWidgets.QGraphicsPathItem):
        def __init__(
            self,
            node_a,
            node_b,
            kind="backbone",
            layer=0,
            lw="",
            linear_layout=False,
        ):
            super(Dssr2DEdgeItem, self).__init__()
            self.node_a = node_a
            self.node_b = node_b
            self.kind = str(kind)
            self.layer = int(layer)
            self.lw = str(lw or "")
            self.linear_layout = bool(linear_layout)
            self.setZValue(-5.0 if self.kind == "backbone" else -3.0)
            self._set_style()
            self.update_geometry()

        def _set_style(self):
            if self.kind == "backbone":
                color = QtGui.QColor(95, 95, 95)
                width = 1.5
                style = QtCore.Qt.SolidLine
            elif self.kind == "tertiary":
                color = QtGui.QColor(145, 90, 170)
                width = 1.25
                style = QtCore.Qt.DashLine
            elif self.layer > 0:
                palette = [
                    QtGui.QColor(210, 60, 60),
                    QtGui.QColor(230, 125, 35),
                    QtGui.QColor(165, 70, 190),
                    QtGui.QColor(30, 150, 150),
                ]
                color = palette[(self.layer - 1) % len(palette)]
                width = 1.8
                style = QtCore.Qt.SolidLine
            else:
                color = QtGui.QColor(60, 95, 205)
                width = 1.8
                style = QtCore.Qt.SolidLine
            pen = QtGui.QPen(color)
            pen.setWidthF(width)
            try:
                pen.setStyle(style)
                pen.setCapStyle(QtCore.Qt.RoundCap)
            except Exception:
                pass
            self.setPen(pen)
            if self.lw:
                self.setToolTip("Base pair: %s" % self.lw)

        def update_geometry(self):
            p1 = self.node_a.pos()
            p2 = self.node_b.pos()
            path = QtGui.QPainterPath()
            path.moveTo(p1)

            if self.linear_layout and self.kind != "backbone":
                span = abs(int(self.node_a.nt_index) - int(self.node_b.nt_index))
                height = min(330.0, 30.0 + 7.0 * span)
                sign = -1.0 if self.layer % 2 else 1.0
                if self.kind == "tertiary":
                    sign *= -1.0
                ctrl = QtCore.QPointF((p1.x() + p2.x()) * 0.5, sign * height)
                path.quadTo(ctrl, p2)
            else:
                path.lineTo(p2)

            self.setPath(path)


    class Dssr2DNodeItem(QtWidgets.QGraphicsEllipseItem):
        RADIUS = 13.0

        def __init__(self, viewer, nt, x, y):
            r = self.RADIUS
            super(Dssr2DNodeItem, self).__init__(-r, -r, 2.0 * r, 2.0 * r)
            self.viewer = viewer
            self.nt = nt
            self.nt_index = int(nt.get("index", 0))
            self.edge_items = []
            self.setPos(float(x), float(y))
            self.setZValue(5.0)

            try:
                flags = (
                    QtWidgets.QGraphicsItem.ItemIsSelectable
                    | QtWidgets.QGraphicsItem.ItemIsMovable
                    | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges
                )
                self.setFlags(flags)
            except Exception:
                pass

            self._apply_style(False)
            self._add_text()
            self.setToolTip(self._tooltip())

        def _base_fill(self):
            base = str(self.nt.get("base", "N")).upper()[:1]
            if not getattr(self.viewer, "base_colors", True):
                return QtGui.QColor(250, 250, 250)
            colors = {
                "A": QtGui.QColor(224, 244, 218),
                "C": QtGui.QColor(218, 235, 251),
                "G": QtGui.QColor(252, 240, 194),
                "U": QtGui.QColor(250, 220, 220),
                "T": QtGui.QColor(250, 220, 220),
                "I": QtGui.QColor(235, 224, 248),
            }
            return colors.get(base, QtGui.QColor(242, 242, 242))

        def _apply_style(self, selected):
            if selected:
                pen = QtGui.QPen(QtGui.QColor(230, 55, 35))
                pen.setWidthF(2.8)
            else:
                pen = QtGui.QPen(QtGui.QColor(45, 45, 45))
                pen.setWidthF(1.2)
            self.setPen(pen)
            self.setBrush(QtGui.QBrush(self._base_fill()))

        def _add_text(self):
            text = QtWidgets.QGraphicsSimpleTextItem(
                str(self.nt.get("base", "N")), self
            )
            font = QtGui.QFont("Sans Serif")
            font.setPointSize(9)
            font.setBold(True)
            text.setFont(font)
            rect = text.boundingRect()
            text.setPos(-rect.width() / 2.0, -rect.height() / 2.0)
            text.setBrush(QtGui.QBrush(QtGui.QColor(25, 25, 25)))
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
            try:
                if change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
                    for edge in list(self.edge_items):
                        edge.update_geometry()
            except Exception:
                pass
            return super(Dssr2DNodeItem, self).itemChange(change, value)

        def mousePressEvent(self, event):
            try:
                self.viewer.select_nucleotide(self.nt_index)
            except Exception:
                pass
            super(Dssr2DNodeItem, self).mousePressEvent(event)

        def set_selected_style(self, selected):
            self._apply_style(bool(selected))


    class DssrPython2DDialog(QtWidgets.QDialog):
        def __init__(
            self,
            model,
            pymol_selection="all",
            algorithm="radiate",
            number_every=10,
            show_tertiary=True,
            parent=None,
        ):
            super(DssrPython2DDialog, self).__init__(parent)
            self.model = model
            self.pymol_selection = str(pymol_selection or "all")
            self.algorithm = str(algorithm or "radiate")
            self.number_every = max(0, int(number_every))
            self.show_tertiary = bool(show_tertiary)
            self.base_colors = True
            self.nodes = []
            self.edges = []
            self._selected_node = None

            self.setWindowTitle("RNA 2D viewer — %s" % self.model.title)
            self.resize(1120, 820)

            root = QtWidgets.QVBoxLayout(self)
            controls = QtWidgets.QHBoxLayout()
            root.addLayout(controls)

            controls.addWidget(QtWidgets.QLabel("layout"))
            self.layout_combo = QtWidgets.QComboBox()
            self.layout_combo.addItems(["radiate", "force", "circular", "linear"])
            idx = self.layout_combo.findText(self.algorithm)
            if idx >= 0:
                self.layout_combo.setCurrentIndex(idx)
            controls.addWidget(self.layout_combo)

            controls.addWidget(QtWidgets.QLabel("number every"))
            self.number_spin = QtWidgets.QSpinBox()
            self.number_spin.setMinimum(0)
            self.number_spin.setMaximum(10000)
            self.number_spin.setValue(self.number_every)
            controls.addWidget(self.number_spin)

            self.tertiary_cb = QtWidgets.QCheckBox("DSSR extra pairs")
            self.tertiary_cb.setChecked(self.show_tertiary)
            controls.addWidget(self.tertiary_cb)

            self.base_colors_cb = QtWidgets.QCheckBox("base colors")
            self.base_colors_cb.setChecked(True)
            controls.addWidget(self.base_colors_cb)

            self.redraw_btn = QtWidgets.QPushButton("redraw")
            self.redraw_btn.clicked.connect(self.redraw)
            controls.addWidget(self.redraw_btn)

            self.fit_btn = QtWidgets.QPushButton("fit")
            self.fit_btn.clicked.connect(self.fit_scene)
            controls.addWidget(self.fit_btn)

            self.copy_btn = QtWidgets.QPushButton("copy DBN")
            self.copy_btn.clicked.connect(self.copy_dbn)
            controls.addWidget(self.copy_btn)

            self.export_btn = QtWidgets.QPushButton("export image")
            self.export_btn.clicked.connect(self.export_image)
            controls.addWidget(self.export_btn)
            controls.addStretch(1)

            self.scene = QtWidgets.QGraphicsScene(self)
            self.view = Dssr2DGraphicsView(self.scene, self)
            root.addWidget(self.view, 1)

            self.status_label = QtWidgets.QLabel(self.model.summary())
            self.status_label.setWordWrap(True)
            root.addWidget(self.status_label)

            self.layout_combo.currentTextChanged.connect(self.redraw)
            self.number_spin.valueChanged.connect(self.redraw)
            self.tertiary_cb.toggled.connect(self.redraw)
            self.base_colors_cb.toggled.connect(self.redraw)

            self.redraw()
            try:
                QtCore.QTimer.singleShot(0, self.fit_scene)
            except Exception:
                self.fit_scene()

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
            return edge

        def redraw(self, *_args):
            algorithm = self.layout_combo.currentText().strip().lower()
            self.algorithm = algorithm
            self.number_every = int(self.number_spin.value())
            self.show_tertiary = bool(self.tertiary_cb.isChecked())
            self.base_colors = bool(self.base_colors_cb.isChecked())

            self.scene.clear()
            self.nodes = []
            self.edges = []
            self._selected_node = None

            positions = Dssr2DLayout.compute(self.model, algorithm)
            linear = algorithm == "linear"

            # Nodes are created before edges so all edges can follow node motion.
            for nt, (x, y) in zip(self.model.nts, positions):
                node = Dssr2DNodeItem(self, nt, x, y)
                self.scene.addItem(node)
                self.nodes.append(node)

            for i in range(len(self.nodes) - 1):
                if i not in self.model.chain_breaks:
                    self._add_edge(i, i + 1, "backbone", linear=linear)

            for pair in self.model.secondary_pairs:
                self._add_edge(
                    int(pair["i"]),
                    int(pair["j"]),
                    "secondary",
                    layer=int(pair.get("layer", 0)),
                    lw=pair.get("lw", ""),
                    linear=linear,
                )

            if self.show_tertiary:
                for pair in self.model.tertiary_pairs:
                    self._add_edge(
                        int(pair["i"]),
                        int(pair["j"]),
                        "tertiary",
                        layer=-1,
                        lw=pair.get("lw", ""),
                        linear=linear,
                    )

            self._add_number_labels()
            self._add_chain_labels()
            rect = self.scene.itemsBoundingRect().adjusted(-70.0, -70.0, 70.0, 70.0)
            self.scene.setSceneRect(rect)
            self.status_label.setText(
                self.model.summary()
                + " | layout=%s | wheel=zoom, drag background=pan, drag bases=edit"
                % algorithm
            )
            self.fit_scene()

        def _add_number_labels(self):
            period = int(self.number_every)
            n = len(self.nodes)
            if n <= 0:
                return
            label_indices = {0, n - 1}
            if period > 0:
                label_indices.update(i for i in range(n) if (i + 1) % period == 0)
            for break_after in self.model.chain_breaks:
                if 0 <= break_after < n:
                    label_indices.add(break_after)
                if 0 <= break_after + 1 < n:
                    label_indices.add(break_after + 1)

            for i in sorted(label_indices):
                node = self.nodes[i]
                nt = self.model.nts[i]
                text_value = str(nt.get("resi") or nt.get("number", i + 1))
                text = QtWidgets.QGraphicsSimpleTextItem(text_value, node)
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(8)
                text.setFont(font)
                text.setBrush(QtGui.QBrush(QtGui.QColor(35, 35, 35)))
                text.setPos(15.0, -24.0)
                text.setZValue(8.0)

        def _add_chain_labels(self):
            seen = set()
            for i, nt in enumerate(self.model.nts):
                chain = str(nt.get("chain", ""))
                if not chain or chain in seen or i >= len(self.nodes):
                    continue
                seen.add(chain)
                label = QtWidgets.QGraphicsSimpleTextItem(
                    "chain %s" % chain, self.nodes[i]
                )
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(9)
                font.setBold(True)
                label.setFont(font)
                label.setBrush(QtGui.QBrush(QtGui.QColor(25, 90, 120)))
                label.setPos(-18.0, -48.0)
                label.setZValue(8.0)

        def fit_scene(self):
            try:
                rect = self.scene.itemsBoundingRect().adjusted(-35, -35, 35, 35)
                self.view.fitInView(rect, QtCore.Qt.KeepAspectRatio)
            except Exception:
                pass

        def _selection_for_nt(self, nt):
            clauses = ["(%s)" % self.pymol_selection]
            chain = str(nt.get("chain", "")).strip()
            resi = str(nt.get("resi", "")).strip()
            if chain:
                clauses.append("chain %s" % chain)
            if resi:
                clauses.append("resi %s" % resi)
            return "byres (%s)" % " and ".join(clauses)

        def select_nucleotide(self, index):
            if index < 0 or index >= len(self.model.nts):
                return
            if self._selected_node is not None:
                try:
                    self._selected_node.set_selected_style(False)
                except Exception:
                    pass
            node = self.nodes[index]
            node.set_selected_style(True)
            self._selected_node = node
            nt = self.model.nts[index]
            expr = self._selection_for_nt(nt)
            try:
                cmd.select("sele", expr)
                self.status_label.setText(
                    "Selected nt %d: %s | %s"
                    % (index + 1, nt.get("nt_id") or expr, self.model.summary())
                )
            except Exception as e:
                self.status_label.setText("PyMOL selection error: %s" % str(e))

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
            image.fill(QtGui.QColor("white"))
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
            rect = self.scene.itemsBoundingRect().adjusted(-30, -30, 30, 30)
            generator = QtSvg.QSvgGenerator()
            generator.setFileName(path)
            generator.setSize(QtCore.QSize(max(1, int(rect.width())), max(1, int(rect.height()))))
            generator.setViewBox(rect)
            generator.setTitle(self.model.title)
            generator.setDescription("RNA secondary structure derived by DSSR")
            painter = QtGui.QPainter(generator)
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            self.scene.render(painter, rect, rect, QtCore.Qt.KeepAspectRatio)
            painter.end()


class DssrPython2DIntegration:
    @staticmethod
    def _open_dialog(
        dssr_data,
        selection="all",
        title="RNA secondary structure",
        algorithm="radiate",
        number_every=10,
        show_tertiary=1,
        parent=None,
    ):
        if QtWidgets is None or QtCore is None or QtGui is None:
            raise CmdException("Qt graphics support is unavailable in this PyMOL build")
        model = Dssr2DModel.from_dssr(dssr_data, title=title)
        dialog = DssrPython2DDialog(
            model=model,
            pymol_selection=selection,
            algorithm=algorithm,
            number_every=number_every,
            show_tertiary=bool(int(show_tertiary)),
            parent=parent,
        )
        _DSSR_PY2D_DIALOGS.append(dialog)

        def _remove_dialog(*_args):
            try:
                _DSSR_PY2D_DIALOGS.remove(dialog)
            except Exception:
                pass

        try:
            dialog.destroyed.connect(_remove_dialog)
        except Exception:
            pass
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        return dialog

    @staticmethod
    def install_gui_controls(dialog):
        if QtWidgets is None or QtCore is None:
            return
        if getattr(dialog, "_dssr_python2d_controls_installed", False):
            return
        dialog._dssr_python2d_controls_installed = True

        group = QtWidgets.QGroupBox("RNA 2D (pure Python — no Java/Jmol)")
        row = QtWidgets.QHBoxLayout(group)

        row.addWidget(QtWidgets.QLabel("layout"))
        dialog.py2d_layout_combo = QtWidgets.QComboBox()
        dialog.py2d_layout_combo.addItems(["radiate", "force", "circular", "linear"])
        row.addWidget(dialog.py2d_layout_combo)

        row.addWidget(QtWidgets.QLabel("number every"))
        dialog.py2d_number_spin = QtWidgets.QSpinBox()
        dialog.py2d_number_spin.setMinimum(0)
        dialog.py2d_number_spin.setMaximum(10000)
        dialog.py2d_number_spin.setValue(10)
        row.addWidget(dialog.py2d_number_spin)

        dialog.py2d_tertiary_cb = QtWidgets.QCheckBox("show DSSR extra pairs")
        dialog.py2d_tertiary_cb.setChecked(True)
        row.addWidget(dialog.py2d_tertiary_cb)

        dialog.py2d_open_btn = QtWidgets.QPushButton("open 2D")
        dialog.py2d_open_btn.setToolTip(
            "Open an in-process Python/Qt RNA secondary-structure viewer. "
            "No Jmol ZIP, Java, browser, or network connection is required."
        )
        dialog.py2d_open_btn.clicked.connect(
            lambda: DssrPython2DIntegration.gui_open(dialog)
        )
        row.addWidget(dialog.py2d_open_btn)
        row.addStretch(1)

        root = dialog.layout()
        if root is not None:
            try:
                root.insertWidget(max(0, root.count() - 1), group)
            except Exception:
                root.addWidget(group)
        dialog.py2d_group = group

    @staticmethod
    def gui_open(dialog):
        try:
            if not cmd.get_object_list():
                raise CmdException(
                    "No structure loaded. Please load a PDB/CIF file before opening 2D view"
                )
            selection, exe, state, precolor_on = dialog._get_dssr_context()
            dssr_data = dialog._get_dssr_data(selection, state, exe, precolor_on)
            algorithm = dialog.py2d_layout_combo.currentText().strip().lower()
            number_every = int(dialog.py2d_number_spin.value())
            show_tertiary = 1 if dialog.py2d_tertiary_cb.isChecked() else 0
            title = "%s state %d — secondary structure derived by DSSR" % (
                selection,
                state,
            )
            viewer = DssrPython2DIntegration._open_dialog(
                dssr_data=dssr_data,
                selection=selection,
                title=title,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=dialog,
            )
            try:
                dialog.status_label.setText(
                    "Opened pure-Python RNA 2D viewer: %s" % viewer.model.summary()
                )
            except Exception:
                pass
            try:
                dialog._append_report(
                    "Python 2D opened: %s; layout=%s"
                    % (viewer.model.summary(), algorithm)
                )
            except Exception:
                pass
        except Exception as e:
            msg = "RNA 2D viewer error: %s" % str(e)
            try:
                dialog.status_label.setText(msg)
            except Exception:
                pass
            try:
                dialog._append_report(msg)
            except Exception:
                pass
            try:
                QtWidgets.QMessageBox.critical(dialog, "RNA 2D viewer", msg)
            except Exception:
                pass


def dssr_2d(
    selection="all",
    state=-1,
    exe="x3dna-dssr",
    layout="radiate",
    number_every=10,
    show_tertiary=1,
    title="",
    quiet=1,
):
    """
    DESCRIPTION

        Open an interactive RNA secondary-structure viewer implemented entirely
        in Python/Qt.  It requires no Java, Jmol, VARNA, ZIP extraction, browser,
        or network connection.  X3DNA-DSSR remains the annotation engine.

    USAGE

        dssr_2d [ selection [, state [, exe [, layout [, number_every
            [, show_tertiary [, title [, quiet ]]]]]]]]

    ARGUMENTS

        selection = str: PyMOL atom selection {default: all}

        state = int: object state {-1: current state}

        exe = str: path to x3dna-dssr {default: x3dna-dssr}

        layout = radiate|force|circular|linear {default: radiate}

        number_every = int: residue-number labeling period; 0 disables periodic
        labels {default: 10}

        show_tertiary = 0|1: show DSSR pairs not present in the secondary DBN
        scaffold as dashed lines {default: 1}

        title = str: optional viewer title

        quiet = 0|1: print a summary {default: 1}

    EXAMPLE

        fetch 1ehz, async=0
        dssr_2d 1ehz
    """
    selection = HelperFunctions.unquote(selection)
    exe = HelperFunctions.unquote(exe)
    layout = HelperFunctions.unquote(layout)
    title = HelperFunctions.unquote(title)
    try:
        state = int(state)
    except Exception:
        state = -1
    try:
        number_every = int(number_every)
    except Exception:
        number_every = 10
    try:
        show_tertiary = int(show_tertiary)
    except Exception:
        show_tertiary = 1
    try:
        quiet = int(quiet)
    except Exception:
        quiet = 1

    if state <= 0:
        try:
            state = int(cmd.get_state())
        except Exception:
            state = 1

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmpfilepdb = tmp.name
    tmp.close()
    try:
        cmd.save(tmpfilepdb, selection, state)
        dssr_data = HelperFunctions.run_dssr_json(tmpfilepdb, exe)
    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass

    if not str(title).strip():
        title = "%s state %d — secondary structure derived by DSSR" % (
            selection,
            state,
        )

    dialog = DssrPython2DIntegration._open_dialog(
        dssr_data=dssr_data,
        selection=selection,
        title=title,
        algorithm=layout,
        number_every=number_every,
        show_tertiary=show_tertiary,
        parent=None,
    )

    if not quiet:
        print("dssr_2d: opened pure-Python viewer — %s" % dialog.model.summary())
    return dialog


# Add the controls without altering the original GUI implementation.
if not getattr(DssrGuiDialog, "_dssr_python2d_patched", False):
    _DSSR_PY2D_ORIGINAL_DIALOG_INIT = DssrGuiDialog.__init__

    def _dssr_python2d_extended_dialog_init(self, *args, **kwargs):
        _DSSR_PY2D_ORIGINAL_DIALOG_INIT(self, *args, **kwargs)
        try:
            DssrPython2DIntegration.install_gui_controls(self)
        except Exception as e:
            try:
                print("DSSR pure-Python 2D GUI setup error: %s" % str(e))
            except Exception:
                pass

    DssrGuiDialog.__init__ = _dssr_python2d_extended_dialog_init
    DssrGuiDialog._dssr_python2d_patched = True


DssrFunctions.dssr_2d = staticmethod(dssr_2d)
cmd.extend("dssr_2d", dssr_2d)

try:
    cmd.auto_arg[0].update({"dssr_2d": cmd.auto_arg[0]["zoom"]})
except Exception:
    pass

try:
    print(
        "Loaded DSSR pure-Python RNA 2D extension %s (no Java/Jmol required)"
        % __DSSR_PY2D_EXTENSION_VERSION__
    )
except Exception:
    pass

# ============================================================================
# DSSR pure-Python RNA 2D publication-layout + manual-editor patch v3.0
#
# This section is APPENDED.  No source code above is deleted or rewritten.
# It adds a topology-aware tRNA cloverleaf layout and a true manual editor:
# free base dragging, group/rubber-band selection, undo/redo, keyboard nudging,
# middle-button or Space+drag panning, and JSON layout save/load.
#
# Design references were studied from VARNA, RNAcanvas, RiboSketch, and Qt's
# Graphics View framework.  This implementation is original Python code and
# does not translate/copy VARNA GPL method bodies.
# ============================================================================

__DSSR_PY2D_EDITOR_PATCH_VERSION__ = "v3.0.0-smart-free-editor"

import bisect


# ---------------------------------------------------------------------------
# Publication-oriented automatic layout
# ---------------------------------------------------------------------------

_DSSR_PY2D_RADIATE_BEFORE_V3 = Dssr2DLayout.radiate
_DSSR_PY2D_COMPUTE_BEFORE_V3 = Dssr2DLayout.compute


def _dssr2d_v3_solve_circle(edge_lengths):
    """Solve a circle whose successive chord lengths close one revolution."""
    lengths = [max(1.0, float(value)) for value in edge_lengths]
    if not lengths:
        return 30.0, []

    lower = max(lengths) * 0.5 + 1.0e-7

    def angles_at(radius):
        return [
            2.0 * math.asin(min(1.0, length / (2.0 * radius)))
            for length in lengths
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


def _dssr2d_v3_normalize(vector, fallback=(0.0, 1.0)):
    length = math.hypot(float(vector[0]), float(vector[1]))
    if length <= 1.0e-12:
        return fallback
    return (float(vector[0]) / length, float(vector[1]) / length)


def _dssr2d_v3_quadratic_equal(indices, start, stop, control, positions):
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
                u * u * start[0]
                + 2.0 * u * t * control[0]
                + t * t * stop[0],
                u * u * start[1]
                + 2.0 * u * t * control[1]
                + t * t * stop[1],
            )
        )

    cumulative = [0.0]
    for first, second in zip(points, points[1:]):
        cumulative.append(
            cumulative[-1]
            + math.hypot(second[0] - first[0], second[1] - first[1])
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


def _dssr2d_v3_place_loop_circle(
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

    outward = _dssr2d_v3_normalize(outward_direction)
    radius, angles = _dssr2d_v3_solve_circle(
        [float(backbone_distance)] * (len(indices) + 1)
        + [float(pair_distance)]
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

    # Both signs close the same polygon; retain the one whose loop bases project
    # farther along the requested outward direction.
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


def _dssr2d_v3_detect_trna_topology(model):
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
        # Each arm needs a non-empty loop or at least enough enclosed residues
        # to behave as a tRNA arm rather than three adjacent bare stems.
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


def _dssr2d_v3_trna_cloverleaf(model, topology):
    """Generate a clean, conventional four-arm tRNA cloverleaf."""
    total = len(model.nts)
    root = topology["root"]
    arms = topology["arms"]
    positions = [None] * total

    pair_distance = 42.0
    helix_rise = 38.0
    backbone_distance = 34.0

    def place_stem(stem, outer_center, direction):
        direction = _dssr2d_v3_normalize(direction)
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

    # Acceptor stem points upward, anticodon arm downward, and the D/T arms
    # spread left/right.  Coordinates are derived from stem lengths and topology.
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
        _dssr2d_v3_place_loop_circle(
            loop_indices,
            positions[arm["inner_i"]],
            positions[arm["inner_j"]],
            direction,
            positions,
            backbone_distance=backbone_distance,
            pair_distance=pair_distance,
        )

    # Four junction/linker arcs connect the acceptor, D, anticodon, and T arms.
    junctions = [
        (root["inner_i"], arms[0]["outer_i"], (-106.0, -106.0)),
        (arms[0]["outer_j"], arms[1]["outer_i"], (-122.0, 96.0)),
        (arms[1]["outer_j"], arms[2]["outer_i"], (122.0, 96.0)),
        (arms[2]["outer_j"], root["inner_j"], (106.0, -106.0)),
    ]
    for first, last, control in junctions:
        _dssr2d_v3_quadratic_equal(
            range(first + 1, last),
            positions[first],
            positions[last],
            control,
            positions,
        )

    # Exterior 5' and 3' tails.  The common 3'-CCA extension is drawn cleanly
    # away from the acceptor stem without assuming any particular sequence.
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

    # Defensive recovery if a modified/atypical tRNA contains a residue outside
    # the four recognized regions.
    fallback = _DSSR_PY2D_RADIATE_BEFORE_V3(model)
    for index in range(total):
        if positions[index] is None:
            positions[index] = fallback[index]

    center_x = sum(point[0] for point in positions) / float(total)
    center_y = sum(point[1] for point in positions) / float(total)
    model._dssr2d_layout_variant = "tRNA cloverleaf"
    return [
        (point[0] - center_x, point[1] - center_y)
        for point in positions
    ]


def _dssr2d_v3_smart_layout(model):
    topology = _dssr2d_v3_detect_trna_topology(model)
    if topology is not None:
        return _dssr2d_v3_trna_cloverleaf(model, topology)
    model._dssr2d_layout_variant = "general radiate"
    return _DSSR_PY2D_RADIATE_BEFORE_V3(model)


def _dssr2d_v3_compute(model, algorithm):
    name = str(algorithm or "smart").strip().lower()
    if name in ("smart", "publication", "auto", "editor"):
        return _dssr2d_v3_smart_layout(model)
    if name in ("radiate", "radial"):
        # Compatibility: radiate also gets the topology-aware tRNA template.
        return _dssr2d_v3_smart_layout(model)
    return _DSSR_PY2D_COMPUTE_BEFORE_V3(model, name)


Dssr2DLayout.smart = staticmethod(_dssr2d_v3_smart_layout)
Dssr2DLayout.radiate = staticmethod(_dssr2d_v3_smart_layout)
Dssr2DLayout.compute = staticmethod(_dssr2d_v3_compute)
Dssr2DLayout._smart_editor_v3_installed = True


# ---------------------------------------------------------------------------
# Curved non-scaffold interaction overlays
# ---------------------------------------------------------------------------

if "Dssr2DEdgeItem" in globals():
    def _dssr2d_v3_edge_geometry(self):
        first = self.node_a.pos()
        second = self.node_b.pos()
        path = QtGui.QPainterPath()
        path.moveTo(first)

        if self.kind == "backbone":
            path.lineTo(second)
        elif self.linear_layout:
            span = abs(int(self.node_a.nt_index) - int(self.node_b.nt_index))
            height = min(350.0, 30.0 + 7.0 * span)
            sign = -1.0 if int(self.layer) % 2 else 1.0
            if self.kind == "tertiary":
                sign *= -1.0
            control = QtCore.QPointF(
                0.5 * (first.x() + second.x()),
                sign * height,
            )
            path.quadTo(control, second)
        elif self.kind == "tertiary" or int(self.layer) > 0:
            dx = second.x() - first.x()
            dy = second.y() - first.y()
            distance = math.hypot(dx, dy)
            if distance <= 1.0e-8:
                path.lineTo(second)
            else:
                nx = -dy / distance
                ny = dx / distance
                middle_x = 0.5 * (first.x() + second.x())
                middle_y = 0.5 * (first.y() + second.y())
                curvature = min(
                    240.0,
                    max(36.0, distance * (0.28 if self.kind == "tertiary" else 0.21)),
                )
                option_a = (
                    middle_x + nx * curvature,
                    middle_y + ny * curvature,
                )
                option_b = (
                    middle_x - nx * curvature,
                    middle_y - ny * curvature,
                )
                if math.hypot(*option_b) > math.hypot(*option_a):
                    nx = -nx
                    ny = -ny
                separation = (
                    int(self.node_a.nt_index) * 17
                    + int(self.node_b.nt_index) * 31
                    + int(self.layer) * 7
                ) % 3
                curvature += (separation - 1) * 12.0
                control = QtCore.QPointF(
                    middle_x + nx * curvature,
                    middle_y + ny * curvature,
                )
                path.quadTo(control, second)
        else:
            path.lineTo(second)
        self.setPath(path)

    Dssr2DEdgeItem.update_geometry = _dssr2d_v3_edge_geometry
    Dssr2DEdgeItem._smart_editor_v3_installed = True


# ---------------------------------------------------------------------------
# True manual editor (Qt Graphics View)
# ---------------------------------------------------------------------------

if QtWidgets is not None and QtCore is not None and QtGui is not None:
    _DSSR_PY2D_VIEW_BEFORE_V3 = Dssr2DGraphicsView
    _DSSR_PY2D_NODE_BEFORE_V3 = Dssr2DNodeItem
    _DSSR_PY2D_DIALOG_BEFORE_V3 = DssrPython2DDialog

    def _dssr2d_v3_no_mouse(item):
        try:
            item.setAcceptedMouseButtons(QtCore.Qt.NoButton)
        except Exception:
            pass

    class Dssr2DGraphicsViewV3(_DSSR_PY2D_VIEW_BEFORE_V3):
        """Rubber-band selection plus dedicated pan/zoom/editor shortcuts."""

        def __init__(self, scene, parent=None):
            super(Dssr2DGraphicsViewV3, self).__init__(scene, parent)
            self.editor = parent
            self._v3_panning = False
            self._v3_space_down = False
            self._v3_pan_last = None
            try:
                self.setInteractive(True)
                self.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
                self.setRubberBandSelectionMode(QtCore.Qt.IntersectsItemShape)
                self.setFocusPolicy(QtCore.Qt.StrongFocus)
                self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
                self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorViewCenter)
            except Exception:
                pass

        def mousePressEvent(self, event):
            try:
                is_middle = event.button() == QtCore.Qt.MiddleButton
                is_space_left = (
                    event.button() == QtCore.Qt.LeftButton and self._v3_space_down
                )
            except Exception:
                is_middle = False
                is_space_left = False

            if is_middle or is_space_left:
                self._v3_panning = True
                self._v3_pan_last = event.pos()
                try:
                    self.setCursor(QtCore.Qt.ClosedHandCursor)
                except Exception:
                    pass
                event.accept()
                return
            super(Dssr2DGraphicsViewV3, self).mousePressEvent(event)

        def mouseMoveEvent(self, event):
            if self._v3_panning and self._v3_pan_last is not None:
                delta = event.pos() - self._v3_pan_last
                self._v3_pan_last = event.pos()
                try:
                    self.horizontalScrollBar().setValue(
                        self.horizontalScrollBar().value() - delta.x()
                    )
                    self.verticalScrollBar().setValue(
                        self.verticalScrollBar().value() - delta.y()
                    )
                except Exception:
                    pass
                event.accept()
                return
            super(Dssr2DGraphicsViewV3, self).mouseMoveEvent(event)

        def mouseReleaseEvent(self, event):
            if self._v3_panning:
                self._v3_panning = False
                self._v3_pan_last = None
                try:
                    self.unsetCursor()
                except Exception:
                    pass
                event.accept()
                return
            super(Dssr2DGraphicsViewV3, self).mouseReleaseEvent(event)

        def keyPressEvent(self, event):
            key = event.key()
            modifiers = event.modifiers()
            control = bool(modifiers & QtCore.Qt.ControlModifier)
            shift = bool(modifiers & QtCore.Qt.ShiftModifier)

            if key == QtCore.Qt.Key_Space:
                self._v3_space_down = True
                try:
                    self.setCursor(QtCore.Qt.OpenHandCursor)
                except Exception:
                    pass
                event.accept()
                return
            if control and key == QtCore.Qt.Key_Z:
                if shift:
                    self.editor.redo_layout()
                else:
                    self.editor.undo_layout()
                event.accept()
                return
            if control and key == QtCore.Qt.Key_Y:
                self.editor.redo_layout()
                event.accept()
                return
            if control and key == QtCore.Qt.Key_A:
                self.editor.select_all_bases()
                event.accept()
                return
            if key == QtCore.Qt.Key_Escape:
                self.editor.clear_base_selection()
                event.accept()
                return
            if key in (
                QtCore.Qt.Key_Left,
                QtCore.Qt.Key_Right,
                QtCore.Qt.Key_Up,
                QtCore.Qt.Key_Down,
            ):
                amount = 10.0 if shift else 2.0
                dx = 0.0
                dy = 0.0
                if key == QtCore.Qt.Key_Left:
                    dx = -amount
                elif key == QtCore.Qt.Key_Right:
                    dx = amount
                elif key == QtCore.Qt.Key_Up:
                    dy = -amount
                elif key == QtCore.Qt.Key_Down:
                    dy = amount
                self.editor.nudge_selected(dx, dy)
                event.accept()
                return
            super(Dssr2DGraphicsViewV3, self).keyPressEvent(event)

        def keyReleaseEvent(self, event):
            if event.key() == QtCore.Qt.Key_Space:
                self._v3_space_down = False
                try:
                    self.unsetCursor()
                except Exception:
                    pass
                event.accept()
                return
            super(Dssr2DGraphicsViewV3, self).keyReleaseEvent(event)

        def mouseDoubleClickEvent(self, event):
            try:
                item = self.itemAt(event.pos())
            except Exception:
                item = None
            if item is None:
                try:
                    self.editor.fit_scene()
                    event.accept()
                    return
                except Exception:
                    pass
            super(Dssr2DGraphicsViewV3, self).mouseDoubleClickEvent(event)


    class Dssr2DNodeItemV3(_DSSR_PY2D_NODE_BEFORE_V3):
        """A nucleotide that can be moved freely, alone or in a selected group."""

        def __init__(self, viewer, nt, x, y):
            super(Dssr2DNodeItemV3, self).__init__(viewer, nt, x, y)
            self._v3_dragging = False
            self._v3_drag_origin = None
            self._v3_drag_starts = {}
            self._v3_drag_before = None
            self._v3_pan_dragging = False
            self._v3_pan_last = None
            try:
                self.setFlag(QtWidgets.QGraphicsItem.ItemIsFocusable, True)
                self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, False)
                self.setAcceptHoverEvents(True)
                self.setCursor(QtCore.Qt.OpenHandCursor)
            except Exception:
                pass
            _dssr2d_v3_no_mouse(getattr(self, "base_text_item", None))

        def itemChange(self, change, value):
            result = super(Dssr2DNodeItemV3, self).itemChange(change, value)
            try:
                if change == QtWidgets.QGraphicsItem.ItemSelectedHasChanged:
                    self._apply_style(bool(value))
                elif change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
                    self.viewer._v3_schedule_scene_rect()
            except Exception:
                pass
            return result

        def mousePressEvent(self, event):
            button = event.button()
            pan_requested = (
                button == QtCore.Qt.MiddleButton
                or (
                    button == QtCore.Qt.LeftButton
                    and bool(getattr(self.viewer.view, "_v3_space_down", False))
                )
            )
            if pan_requested:
                self._v3_pan_dragging = True
                try:
                    self._v3_pan_last = event.screenPos()
                    self.viewer.view.setCursor(QtCore.Qt.ClosedHandCursor)
                except Exception:
                    self._v3_pan_last = None
                event.accept()
                return

            if button != QtCore.Qt.LeftButton:
                super(Dssr2DNodeItemV3, self).mousePressEvent(event)
                return

            modifiers = event.modifiers()
            control = bool(modifiers & QtCore.Qt.ControlModifier)
            shift = bool(modifiers & QtCore.Qt.ShiftModifier)
            scene = self.scene()

            if control:
                self.setSelected(not self.isSelected())
                if not self.isSelected():
                    self.viewer._v3_active_index = None
                    self.viewer._v3_sync_pymol_selection()
                    event.accept()
                    return
            elif shift:
                self.setSelected(True)
            else:
                if not self.isSelected():
                    try:
                        scene.clearSelection()
                    except Exception:
                        pass
                    self.setSelected(True)

            self.viewer._v3_active_index = self.nt_index
            selected = [node for node in self.viewer.nodes if node.isSelected()]
            if self not in selected:
                selected.append(self)
                self.setSelected(True)

            self._v3_dragging = True
            self._v3_drag_origin = event.scenePos()
            self._v3_drag_before = self.viewer._v3_capture_positions()
            self._v3_drag_starts = {
                node.nt_index: QtCore.QPointF(node.pos()) for node in selected
            }
            try:
                self.setCursor(QtCore.Qt.ClosedHandCursor)
                self.viewer.view.setFocus()
            except Exception:
                pass
            self.viewer._v3_sync_pymol_selection()
            event.accept()

        def mouseMoveEvent(self, event):
            if self._v3_pan_dragging and self._v3_pan_last is not None:
                current = event.screenPos()
                delta = current - self._v3_pan_last
                self._v3_pan_last = current
                view = self.viewer.view
                view.horizontalScrollBar().setValue(
                    view.horizontalScrollBar().value() - delta.x()
                )
                view.verticalScrollBar().setValue(
                    view.verticalScrollBar().value() - delta.y()
                )
                event.accept()
                return

            if not self._v3_dragging or self._v3_drag_origin is None:
                super(Dssr2DNodeItemV3, self).mouseMoveEvent(event)
                return

            delta = event.scenePos() - self._v3_drag_origin
            snap = False
            grid = 1.0
            try:
                snap = bool(self.viewer.snap_cb.isChecked())
                grid = max(1.0, float(self.viewer.grid_spin.value()))
            except Exception:
                pass

            for index, start in self._v3_drag_starts.items():
                if index < 0 or index >= len(self.viewer.nodes):
                    continue
                x = start.x() + delta.x()
                y = start.y() + delta.y()
                if snap:
                    x = round(x / grid) * grid
                    y = round(y / grid) * grid
                self.viewer.nodes[index].setPos(x, y)
            event.accept()

        def mouseReleaseEvent(self, event):
            if self._v3_pan_dragging:
                self._v3_pan_dragging = False
                self._v3_pan_last = None
                try:
                    self.viewer.view.unsetCursor()
                except Exception:
                    pass
                event.accept()
                return

            if self._v3_dragging:
                self._v3_dragging = False
                after = self.viewer._v3_capture_positions()
                self.viewer._v3_push_history(
                    self._v3_drag_before,
                    after,
                    "move base%s"
                    % ("s" if len(self._v3_drag_starts) != 1 else ""),
                )
                self._v3_drag_origin = None
                self._v3_drag_starts = {}
                self._v3_drag_before = None
                try:
                    self.setCursor(QtCore.Qt.OpenHandCursor)
                except Exception:
                    pass
                self.viewer._v3_update_editor_status("manual move")
                event.accept()
                return
            super(Dssr2DNodeItemV3, self).mouseReleaseEvent(event)

        def contextMenuEvent(self, event):
            menu = QtWidgets.QMenu()
            select_only = menu.addAction("Select only this base")
            add_pair = menu.addAction("Select its base pair")
            add_stem = menu.addAction("Select its whole stem")
            menu.addSeparator()
            reset_selected = menu.addAction("Reset selected bases to automatic layout")
            undo_action = menu.addAction("Undo last move")

            try:
                chosen = menu.exec_(event.screenPos())
            except Exception:
                try:
                    chosen = menu.exec(event.screenPos())
                except Exception:
                    chosen = None

            if chosen == select_only:
                self.viewer.select_indices([self.nt_index], replace=True)
            elif chosen == add_pair:
                indices = [self.nt_index]
                partner = self.viewer._v3_pair_partner(self.nt_index)
                if partner >= 0:
                    indices.append(partner)
                self.viewer.select_indices(indices, replace=True)
            elif chosen == add_stem:
                self.viewer.select_indices(
                    self.viewer._v3_stem_indices(self.nt_index), replace=True
                )
            elif chosen == reset_selected:
                self.viewer.reset_selected_bases()
            elif chosen == undo_action:
                self.viewer.undo_layout()
            event.accept()


    # Runtime lookups in the original dialog now create the enhanced classes.
    Dssr2DGraphicsView = Dssr2DGraphicsViewV3
    Dssr2DNodeItem = Dssr2DNodeItemV3


    class DssrPython2DDialogV3(_DSSR_PY2D_DIALOG_BEFORE_V3):
        HISTORY_LIMIT = 100

        def __init__(
            self,
            model,
            pymol_selection="all",
            algorithm="smart",
            number_every=10,
            show_tertiary=False,
            parent=None,
        ):
            self._v3_rebuilding = False
            self._v3_force_relayout = True
            self._v3_auto_positions = []
            self._v3_undo = []
            self._v3_redo = []
            self._v3_active_index = None
            self._v3_scene_rect_pending = False
            self._v3_selection_sync_pending = False
            self._v3_editor_installed = False
            if str(algorithm or "").strip().lower() in ("", "radiate"):
                algorithm = "smart"
            self._v3_requested_algorithm = str(algorithm or "smart").strip().lower()

            super(DssrPython2DDialogV3, self).__init__(
                model=model,
                pymol_selection=pymol_selection,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=parent,
            )
            self._v3_install_editor_controls()
            try:
                self.scene.selectionChanged.connect(self._v3_selection_changed)
            except Exception:
                pass
            self._v3_update_history_buttons()
            self._v3_update_editor_status("ready")

        # ---------- scene construction and persistence ----------

        def redraw(self, *_args):
            if self._v3_rebuilding:
                return
            self._v3_rebuilding = True
            try:
                algorithm = self.layout_combo.currentText().strip().lower() or "smart"
                self.algorithm = algorithm
                self.number_every = int(self.number_spin.value())
                self.show_tertiary = bool(self.tertiary_cb.isChecked())
                self.base_colors = bool(self.base_colors_cb.isChecked())

                sender = self.sender()
                current_positions = self._v3_capture_positions()
                selected_indices = [
                    node.nt_index for node in self.nodes if node.isSelected()
                ]
                preserve_positions = (
                    bool(current_positions)
                    and not self._v3_force_relayout
                    and sender is not self.layout_combo
                    and sender is not self.redraw_btn
                )

                old_transform = None
                old_center = None
                if preserve_positions:
                    try:
                        old_transform = QtGui.QTransform(self.view.transform())
                        old_center = self.view.mapToScene(
                            self.view.viewport().rect().center()
                        )
                    except Exception:
                        old_transform = None
                        old_center = None

                if preserve_positions:
                    positions = current_positions
                else:
                    positions = Dssr2DLayout.compute(self.model, algorithm)
                    self._v3_auto_positions = [
                        (float(x), float(y)) for x, y in positions
                    ]

                try:
                    self.scene.blockSignals(True)
                except Exception:
                    pass
                self.scene.clear()
                self.nodes = []
                self.edges = []
                self._selected_node = None

                for nt, (x, y) in zip(self.model.nts, positions):
                    node = Dssr2DNodeItem(self, nt, x, y)
                    self.scene.addItem(node)
                    self.nodes.append(node)

                linear = algorithm == "linear"
                for index in range(len(self.nodes) - 1):
                    if index not in self.model.chain_breaks:
                        self._add_edge(index, index + 1, "backbone", linear=linear)

                for pair in self.model.secondary_pairs:
                    self._add_edge(
                        int(pair["i"]),
                        int(pair["j"]),
                        "secondary",
                        layer=int(pair.get("layer", 0)),
                        lw=pair.get("lw", ""),
                        linear=linear,
                    )

                if self.show_tertiary:
                    for pair in self.model.tertiary_pairs:
                        self._add_edge(
                            int(pair["i"]),
                            int(pair["j"]),
                            "tertiary",
                            layer=-1,
                            lw=pair.get("lw", ""),
                            linear=linear,
                        )

                self._add_number_labels()
                self._add_chain_labels()
                self._v3_update_scene_rect()

                for index in selected_indices:
                    if 0 <= index < len(self.nodes):
                        self.nodes[index].setSelected(True)
                try:
                    self.scene.blockSignals(False)
                except Exception:
                    pass

                if preserve_positions and old_transform is not None:
                    try:
                        self.view.setTransform(old_transform)
                        if old_center is not None:
                            self.view.centerOn(old_center)
                    except Exception:
                        pass
                else:
                    self.fit_scene()

                self._v3_force_relayout = False
                self._v3_update_editor_status(
                    "redrawn" if preserve_positions else "automatic layout"
                )
            finally:
                try:
                    self.scene.blockSignals(False)
                except Exception:
                    pass
                self._v3_rebuilding = False

        def _add_number_labels(self):
            total = len(self.nodes)
            if total <= 0:
                return
            period = int(self.number_every)
            indices = {0, total - 1}
            if period > 0:
                indices.update(
                    index for index in range(total) if (index + 1) % period == 0
                )
            for break_after in self.model.chain_breaks:
                if 0 <= break_after < total:
                    indices.add(break_after)
                if 0 <= break_after + 1 < total:
                    indices.add(break_after + 1)

            center_x = sum(node.pos().x() for node in self.nodes) / float(total)
            center_y = sum(node.pos().y() for node in self.nodes) / float(total)
            for index in sorted(indices):
                node = self.nodes[index]
                nt = self.model.nts[index]
                value = str(nt.get("resi") or nt.get("number", index + 1))
                label = QtWidgets.QGraphicsSimpleTextItem(value, node)
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(8)
                label.setFont(font)
                label.setBrush(QtGui.QBrush(QtGui.QColor(35, 35, 35)))
                dx = node.pos().x() - center_x
                dy = node.pos().y() - center_y
                length = math.hypot(dx, dy)
                if length <= 1.0e-8:
                    dx, dy, length = 1.0, -1.0, math.sqrt(2.0)
                dx = 25.0 * dx / length
                dy = 25.0 * dy / length
                rect = label.boundingRect()
                label.setPos(dx - 0.5 * rect.width(), dy - 0.5 * rect.height())
                label.setZValue(8.0)
                _dssr2d_v3_no_mouse(label)

        def _add_chain_labels(self):
            total = len(self.nodes)
            if total <= 0:
                return
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

            for segment_number, (first, last) in enumerate(segments, 1):
                for index, text_value, offset in (
                    (first, "5′", (-35.0, -22.0)),
                    (last, "3′", (20.0, -22.0)),
                ):
                    label = QtWidgets.QGraphicsSimpleTextItem(
                        text_value, self.nodes[index]
                    )
                    font = QtGui.QFont("Sans Serif")
                    font.setPointSize(10)
                    font.setBold(True)
                    label.setFont(font)
                    label.setBrush(QtGui.QBrush(QtGui.QColor(25, 90, 120)))
                    label.setPos(offset[0], offset[1])
                    label.setZValue(8.0)
                    _dssr2d_v3_no_mouse(label)

                if len(segments) > 1:
                    chain = str(self.model.nts[first].get("chain", "")).strip()
                    text_value = "chain %s" % (chain or segment_number)
                    label = QtWidgets.QGraphicsSimpleTextItem(
                        text_value, self.nodes[first]
                    )
                    font = QtGui.QFont("Sans Serif")
                    font.setPointSize(8)
                    font.setBold(True)
                    label.setFont(font)
                    label.setBrush(QtGui.QBrush(QtGui.QColor(25, 90, 120)))
                    label.setPos(-42.0, -46.0)
                    label.setZValue(8.0)
                    _dssr2d_v3_no_mouse(label)

        def _v3_capture_positions(self):
            if not getattr(self, "nodes", None):
                return []
            return [
                (float(node.pos().x()), float(node.pos().y()))
                for node in self.nodes
            ]

        def _v3_apply_positions(self, positions):
            if len(positions) != len(self.nodes):
                raise CmdException(
                    "Layout has %d coordinates but this RNA has %d nucleotides"
                    % (len(positions), len(self.nodes))
                )
            for node, point in zip(self.nodes, positions):
                node.setPos(float(point[0]), float(point[1]))
            for edge in self.edges:
                try:
                    edge.update_geometry()
                except Exception:
                    pass
            self._v3_update_scene_rect()

        def _v3_positions_changed(self, first, second):
            if not first or not second or len(first) != len(second):
                return False
            return any(
                abs(a[0] - b[0]) > 1.0e-5 or abs(a[1] - b[1]) > 1.0e-5
                for a, b in zip(first, second)
            )

        def _v3_push_history(self, before, after, label="edit"):
            if not self._v3_positions_changed(before, after):
                return
            self._v3_undo.append(
                {
                    "before": [(float(x), float(y)) for x, y in before],
                    "after": [(float(x), float(y)) for x, y in after],
                    "label": str(label),
                }
            )
            if len(self._v3_undo) > self.HISTORY_LIMIT:
                self._v3_undo = self._v3_undo[-self.HISTORY_LIMIT :]
            self._v3_redo = []
            self._v3_update_history_buttons()

        def undo_layout(self):
            if not self._v3_undo:
                return
            entry = self._v3_undo.pop()
            self._v3_apply_positions(entry["before"])
            self._v3_redo.append(entry)
            self._v3_update_history_buttons()
            self._v3_update_editor_status("undo: %s" % entry["label"])

        def redo_layout(self):
            if not self._v3_redo:
                return
            entry = self._v3_redo.pop()
            self._v3_apply_positions(entry["after"])
            self._v3_undo.append(entry)
            self._v3_update_history_buttons()
            self._v3_update_editor_status("redo: %s" % entry["label"])

        def reset_layout(self):
            before = self._v3_capture_positions()
            algorithm = self.layout_combo.currentText().strip().lower() or "smart"
            after = Dssr2DLayout.compute(self.model, algorithm)
            self._v3_auto_positions = [(float(x), float(y)) for x, y in after]
            self._v3_apply_positions(after)
            self._v3_push_history(before, self._v3_capture_positions(), "reset layout")
            self.fit_scene()
            self._v3_update_editor_status("automatic layout reset")

        def reset_selected_bases(self):
            if len(self._v3_auto_positions) != len(self.nodes):
                self._v3_auto_positions = Dssr2DLayout.compute(
                    self.model,
                    self.layout_combo.currentText().strip().lower() or "smart",
                )
            selected = [node for node in self.nodes if node.isSelected()]
            if not selected:
                return
            before = self._v3_capture_positions()
            for node in selected:
                x, y = self._v3_auto_positions[node.nt_index]
                node.setPos(float(x), float(y))
            after = self._v3_capture_positions()
            self._v3_push_history(before, after, "reset selected")
            self._v3_update_scene_rect()
            self._v3_update_editor_status("selected bases reset")

        def nudge_selected(self, dx, dy):
            selected = [node for node in self.nodes if node.isSelected()]
            if not selected:
                return
            before = self._v3_capture_positions()
            for node in selected:
                node.moveBy(float(dx), float(dy))
            after = self._v3_capture_positions()
            self._v3_push_history(before, after, "nudge")
            self._v3_update_editor_status("nudged %d base(s)" % len(selected))

        # ---------- selection ----------

        def select_indices(self, indices, replace=True):
            wanted = set(int(index) for index in indices)
            if replace:
                try:
                    self.scene.clearSelection()
                except Exception:
                    pass
            for index in wanted:
                if 0 <= index < len(self.nodes):
                    self.nodes[index].setSelected(True)
            self._v3_active_index = min(wanted) if wanted else None
            self._v3_sync_pymol_selection()

        def select_all_bases(self):
            self.select_indices(range(len(self.nodes)), replace=True)

        def clear_base_selection(self):
            try:
                self.scene.clearSelection()
            except Exception:
                pass
            self._v3_active_index = None
            self._v3_sync_pymol_selection()

        def select_nucleotide(self, index):
            self.select_indices([index], replace=True)

        def _v3_pair_partner(self, index):
            table = Dssr2DLayout._planar_pair_table(self.model)
            if 0 <= index < len(table):
                return int(table[index])
            return -1

        def _v3_stem_indices(self, index):
            table = Dssr2DLayout._planar_pair_table(self.model)
            stems = Dssr2DLayout._stem_tree(table, self.model.chain_breaks)
            for stem in stems:
                indices = set()
                for first, second in stem.get("pairs", []):
                    indices.add(int(first))
                    indices.add(int(second))
                if index in indices:
                    return sorted(indices)
            partner = self._v3_pair_partner(index)
            return [index] if partner < 0 else sorted([index, partner])

        def _v3_selection_changed(self):
            if self._v3_rebuilding:
                return
            self._v3_sync_pymol_selection()
            self._v3_update_editor_status("selection changed")

        def _v3_sync_pymol_selection(self):
            selected = sorted(
                (node for node in self.nodes if node.isSelected()),
                key=lambda node: node.nt_index,
            )
            if not selected:
                try:
                    cmd.select("sele", "none")
                except Exception:
                    pass
                return

            residues = set()
            fallback_clauses = []
            for node in selected:
                nt = self.model.nts[node.nt_index]
                chain = str(nt.get("chain", "")).strip()
                resi = str(nt.get("resi", "")).strip()
                if chain and resi:
                    residues.add((chain, resi))
                elif resi:
                    fallback_clauses.append("resi %s" % resi)

            core = ParsingAlgos._compact_sel_from_residues(residues)
            pieces = []
            if core:
                pieces.append("(%s)" % core)
            pieces.extend("(%s)" % clause for clause in fallback_clauses)
            if pieces:
                expression = "((%s) and (%s))" % (
                    self.pymol_selection,
                    " or ".join(pieces),
                )
                try:
                    cmd.select("sele", "byres (%s)" % expression)
                except Exception as error:
                    self.status_label.setText(
                        "PyMOL selection error: %s" % str(error)
                    )

        # ---------- layout files ----------

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
                    [round(x, 6), round(y, 6)]
                    for x, y in self._v3_capture_positions()
                ],
            }
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, ensure_ascii=False)
            self._v3_update_editor_status("layout saved: %s" % path)

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
            before = self._v3_capture_positions()
            self._v3_apply_positions(positions)
            self._v3_push_history(before, self._v3_capture_positions(), "load layout")
            self.fit_scene()
            self._v3_update_editor_status("layout loaded: %s" % path)

        # ---------- UI/state helpers ----------

        def _v3_install_editor_controls(self):
            if self._v3_editor_installed:
                return
            self._v3_editor_installed = True

            if self.layout_combo.findText("smart") < 0:
                self.layout_combo.insertItem(0, "smart")
            preferred = str(
                getattr(self, "_v3_requested_algorithm", self.algorithm)
            ).strip().lower()
            index = self.layout_combo.findText(preferred)
            if index < 0:
                index = self.layout_combo.findText("smart")
            if index >= 0:
                self.layout_combo.setCurrentIndex(index)

            try:
                self.redraw_btn.clicked.disconnect()
            except Exception:
                pass
            self.redraw_btn.setText("reset layout")
            self.redraw_btn.setToolTip(
                "Recompute the automatic layout. Manual coordinates are otherwise "
                "preserved when labels/colors/extra-pairs are toggled."
            )
            self.redraw_btn.clicked.connect(self.reset_layout)

            group = QtWidgets.QGroupBox("manual editor")
            row = QtWidgets.QHBoxLayout(group)

            self.undo_btn = QtWidgets.QPushButton("undo")
            self.undo_btn.clicked.connect(self.undo_layout)
            row.addWidget(self.undo_btn)

            self.redo_btn = QtWidgets.QPushButton("redo")
            self.redo_btn.clicked.connect(self.redo_layout)
            row.addWidget(self.redo_btn)

            self.select_all_btn = QtWidgets.QPushButton("select all")
            self.select_all_btn.clicked.connect(self.select_all_bases)
            row.addWidget(self.select_all_btn)

            self.clear_selection_btn = QtWidgets.QPushButton("clear selection")
            self.clear_selection_btn.clicked.connect(self.clear_base_selection)
            row.addWidget(self.clear_selection_btn)

            self.reset_selected_btn = QtWidgets.QPushButton("reset selected")
            self.reset_selected_btn.clicked.connect(self.reset_selected_bases)
            row.addWidget(self.reset_selected_btn)

            self.save_layout_btn = QtWidgets.QPushButton("save layout")
            self.save_layout_btn.clicked.connect(self.save_layout)
            row.addWidget(self.save_layout_btn)

            self.load_layout_btn = QtWidgets.QPushButton("load layout")
            self.load_layout_btn.clicked.connect(self.load_layout)
            row.addWidget(self.load_layout_btn)

            self.snap_cb = QtWidgets.QCheckBox("snap")
            self.snap_cb.setChecked(False)
            self.snap_cb.setToolTip(
                "Optional coordinate grid. Off by default: bases move completely freely."
            )
            row.addWidget(self.snap_cb)

            self.grid_spin = QtWidgets.QDoubleSpinBox()
            self.grid_spin.setMinimum(1.0)
            self.grid_spin.setMaximum(100.0)
            self.grid_spin.setSingleStep(1.0)
            self.grid_spin.setValue(5.0)
            self.grid_spin.setSuffix(" px")
            self.grid_spin.setEnabled(False)
            self.snap_cb.toggled.connect(self.grid_spin.setEnabled)
            row.addWidget(self.grid_spin)

            help_label = QtWidgets.QLabel(
                "drag base • Shift/Ctrl multi-select • drag empty area=box select • "
                "middle or Space+drag=pan • arrows=nudge"
            )
            help_label.setWordWrap(True)
            row.addWidget(help_label, 1)

            root = self.layout()
            if root is not None:
                try:
                    root.insertWidget(1, group)
                except Exception:
                    root.addWidget(group)
            self.editor_group = group

        def _v3_update_history_buttons(self):
            try:
                self.undo_btn.setEnabled(bool(self._v3_undo))
                self.redo_btn.setEnabled(bool(self._v3_redo))
                if self._v3_undo:
                    self.undo_btn.setToolTip(
                        "Undo: %s" % self._v3_undo[-1].get("label", "edit")
                    )
                if self._v3_redo:
                    self.redo_btn.setToolTip(
                        "Redo: %s" % self._v3_redo[-1].get("label", "edit")
                    )
            except Exception:
                pass

        def _v3_update_editor_status(self, action=""):
            selected = sum(1 for node in self.nodes if node.isSelected())
            variant = str(
                getattr(self.model, "_dssr2d_layout_variant", self.algorithm)
            )
            text = (
                self.model.summary()
                + " | layout=%s (%s)" % (self.algorithm, variant)
                + " | selected=%d" % selected
                + " | drag any base freely; Shift/Ctrl=multi-select; middle/Space=pan"
            )
            if action:
                text += " | %s" % action
            self.status_label.setText(text)

        def _v3_schedule_scene_rect(self):
            if self._v3_scene_rect_pending:
                return
            self._v3_scene_rect_pending = True

            def update():
                self._v3_scene_rect_pending = False
                self._v3_update_scene_rect()

            try:
                QtCore.QTimer.singleShot(0, update)
            except Exception:
                update()

        def _v3_update_scene_rect(self):
            try:
                rect = self.scene.itemsBoundingRect().adjusted(
                    -120.0, -120.0, 120.0, 120.0
                )
                self.scene.setSceneRect(rect)
            except Exception:
                pass


    DssrPython2DDialog = DssrPython2DDialogV3
    DssrPython2DDialog._smart_editor_v3_installed = True


# ---------------------------------------------------------------------------
# Main DSSR GUI and PyMOL command defaults
# ---------------------------------------------------------------------------

if not getattr(DssrPython2DIntegration, "_smart_editor_v3_installed", False):
    _DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V3 = (
        DssrPython2DIntegration.install_gui_controls
    )
    _DSSR_PY2D_OPEN_DIALOG_BEFORE_V3 = DssrPython2DIntegration._open_dialog

    def _dssr2d_v3_install_gui_controls(dialog):
        _DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V3(dialog)
        try:
            combo = dialog.py2d_layout_combo
            if combo.findText("smart") < 0:
                combo.insertItem(0, "smart")
            combo.setCurrentIndex(combo.findText("smart"))
            dialog.py2d_tertiary_cb.setChecked(False)
            dialog.py2d_tertiary_cb.setToolTip(
                "Optional non-secondary DSSR contacts. Off by default so the "
                "secondary-structure scaffold stays readable."
            )
            dialog.py2d_open_btn.setToolTip(
                "Open the pure-Python interactive editor. Bases can be dragged "
                "freely; no Java, Jmol, browser, or network is required."
            )
        except Exception:
            pass

    def _dssr2d_v3_open_dialog(
        dssr_data,
        selection="all",
        title="RNA secondary structure",
        algorithm="smart",
        number_every=10,
        show_tertiary=0,
        parent=None,
    ):
        return _DSSR_PY2D_OPEN_DIALOG_BEFORE_V3(
            dssr_data=dssr_data,
            selection=selection,
            title=title,
            algorithm=algorithm,
            number_every=number_every,
            show_tertiary=show_tertiary,
            parent=parent,
        )

    DssrPython2DIntegration.install_gui_controls = staticmethod(
        _dssr2d_v3_install_gui_controls
    )
    DssrPython2DIntegration._open_dialog = staticmethod(
        _dssr2d_v3_open_dialog
    )
    DssrPython2DIntegration._smart_editor_v3_installed = True


_DSSR_PY2D_COMMAND_BEFORE_V3 = dssr_2d


def dssr_2d(
    selection="all",
    state=-1,
    exe="x3dna-dssr",
    layout="smart",
    number_every=10,
    show_tertiary=0,
    title="",
    quiet=1,
):
    """Open the pure-Python RNA 2D publication layout and manual editor.

    ``layout=smart`` recognizes tRNA-like topology and otherwise uses the
    general radiate layout.  Every nucleotide can be dragged independently;
    selected groups move together.  No Java/Jmol dependency is used.
    """
    return _DSSR_PY2D_COMMAND_BEFORE_V3(
        selection=selection,
        state=state,
        exe=exe,
        layout=layout,
        number_every=number_every,
        show_tertiary=show_tertiary,
        title=title,
        quiet=quiet,
    )


DssrFunctions.dssr_2d = staticmethod(dssr_2d)
cmd.extend("dssr_2d", dssr_2d)

try:
    print(
        "Loaded DSSR Python RNA 2D editor %s: smart tRNA layout; free base "
        "dragging; multi-select; undo/redo; save/load"
        % __DSSR_PY2D_EDITOR_PATCH_VERSION__
    )
except Exception:
    pass

# ============================================================================
# DSSR pure-Python RNA 2D editor hotfix v3.0.1
#
# APPEND-ONLY HOTFIX: nothing above is removed or rewritten.
#
# Why this is needed:
# The v3 editor replaces the public class names (DssrPython2DDialog,
# Dssr2DGraphicsView, Dssr2DNodeItem) with enhanced subclasses.  Some legacy
# methods use the two-argument spelling ``super(ClassName, self)``.  Because
# ``ClassName`` is looked up from the module globals at call time, rebinding the
# public name to the subclass can make the legacy method resolve to itself
# again, producing ``RecursionError: maximum recursion depth exceeded`` during
# dialog/node construction.  The wrappers below temporarily bind each legacy
# method's original class name while that method runs.
# ============================================================================

__DSSR_PY2D_EDITOR_HOTFIX_VERSION__ = "v3.0.1-bound-super-hotfix"

import functools as _dssr2d_functools
import traceback as _dssr2d_traceback


def _dssr2d_bind_legacy_super(owner_class, method_name, global_class_name):
    """Make a legacy ``super(ClassName, self)`` method alias-safe.

    This does not alter the method body.  It only restores the class-name
    binding expected by that method for the duration of the call, then restores
    the enhanced public class alias.  PyMOL's GUI runs these calls on the Qt GUI
    thread, so the temporary binding is deterministic in normal plugin use.
    """
    if owner_class is None:
        return False
    try:
        method = getattr(owner_class, method_name)
    except Exception:
        return False
    if getattr(method, "_dssr2d_bound_super_hotfix", False):
        return True

    @_dssr2d_functools.wraps(method)
    def alias_safe_method(*args, **kwargs):
        namespace = globals()
        marker = object()
        previous = namespace.get(global_class_name, marker)
        namespace[global_class_name] = owner_class
        try:
            return method(*args, **kwargs)
        finally:
            if previous is marker:
                namespace.pop(global_class_name, None)
            else:
                namespace[global_class_name] = previous

    alias_safe_method._dssr2d_bound_super_hotfix = True
    alias_safe_method._dssr2d_original_method = method
    setattr(owner_class, method_name, alias_safe_method)
    return True


_DSSR_PY2D_HOTFIXED_METHODS = []

if QtWidgets is not None and QtCore is not None and QtGui is not None:
    # Legacy dialog constructor: fixes super(DssrPython2DDialog, self).
    if "_DSSR_PY2D_DIALOG_BEFORE_V3" in globals():
        if _dssr2d_bind_legacy_super(
            _DSSR_PY2D_DIALOG_BEFORE_V3,
            "__init__",
            "DssrPython2DDialog",
        ):
            _DSSR_PY2D_HOTFIXED_METHODS.append("dialog.__init__")

    # Legacy graphics-view methods: fixes super(Dssr2DGraphicsView, self).
    if "_DSSR_PY2D_VIEW_BEFORE_V3" in globals():
        for _method_name in ("__init__", "wheelEvent"):
            if _dssr2d_bind_legacy_super(
                _DSSR_PY2D_VIEW_BEFORE_V3,
                _method_name,
                "Dssr2DGraphicsView",
            ):
                _DSSR_PY2D_HOTFIXED_METHODS.append(
                    "view.%s" % _method_name
                )

    # Legacy nucleotide methods: fixes super(Dssr2DNodeItem, self).
    if "_DSSR_PY2D_NODE_BEFORE_V3" in globals():
        for _method_name in ("__init__", "itemChange", "mousePressEvent"):
            if _dssr2d_bind_legacy_super(
                _DSSR_PY2D_NODE_BEFORE_V3,
                _method_name,
                "Dssr2DNodeItem",
            ):
                _DSSR_PY2D_HOTFIXED_METHODS.append(
                    "node.%s" % _method_name
                )


# Replace the shallow GUI error message with a complete diagnostic traceback.
# This is intentionally appended as a runtime override; the original method is
# retained above unchanged.
if "DssrPython2DIntegration" in globals():

    @staticmethod
    def _dssr2d_gui_open_with_traceback(dialog):
        try:
            if not cmd.get_object_list():
                raise CmdException(
                    "No structure loaded. Please load a PDB/CIF file before "
                    "opening 2D view"
                )

            selection, exe, state, precolor_on = dialog._get_dssr_context()
            dssr_data = dialog._get_dssr_data(
                selection, state, exe, precolor_on
            )
            algorithm = (
                dialog.py2d_layout_combo.currentText().strip().lower()
                or "smart"
            )
            number_every = int(dialog.py2d_number_spin.value())
            show_tertiary = 1 if dialog.py2d_tertiary_cb.isChecked() else 0
            title = "%s state %d — secondary structure derived by DSSR" % (
                selection,
                state,
            )

            viewer = DssrPython2DIntegration._open_dialog(
                dssr_data=dssr_data,
                selection=selection,
                title=title,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=dialog,
            )
            try:
                dialog.status_label.setText(
                    "Opened pure-Python RNA 2D editor: %s"
                    % viewer.model.summary()
                )
            except Exception:
                pass
            try:
                dialog._append_report(
                    "Python 2D opened: %s; layout=%s; editor hotfix=%s"
                    % (
                        viewer.model.summary(),
                        algorithm,
                        __DSSR_PY2D_EDITOR_HOTFIX_VERSION__,
                    )
                )
            except Exception:
                pass
            return viewer

        except Exception as error:
            trace = _dssr2d_traceback.format_exc()
            msg = "RNA 2D editor error: %s" % str(error)
            try:
                print(msg)
                print(trace)
            except Exception:
                pass
            try:
                dialog.status_label.setText(msg)
            except Exception:
                pass
            try:
                dialog._append_report(msg + "\n\n" + trace)
            except Exception:
                pass
            try:
                QtWidgets.QMessageBox.critical(
                    dialog,
                    "RNA 2D editor",
                    msg
                    + "\n\nA complete traceback was printed in the PyMOL "
                    "console and added to the report panel.",
                )
            except Exception:
                pass
            return None

    DssrPython2DIntegration.gui_open = _dssr2d_gui_open_with_traceback


try:
    print(
        "Loaded DSSR RNA 2D editor hotfix %s (%s)"
        % (
            __DSSR_PY2D_EDITOR_HOTFIX_VERSION__,
            ", ".join(_DSSR_PY2D_HOTFIXED_METHODS) or "no Qt methods patched",
        )
    )
except Exception:
    pass

# ============================================================================
# DSSR pure-Python RNA 2D editor v4.0 — Gel UI + structure-aware interaction
#
# APPEND-ONLY PATCH: all code above remains intact.
#
# Design goals:
#   * glossy, soft "gel" nucleotides with springy hover/press feedback
#   * left-drag on empty canvas pans; Shift/Ctrl-drag makes box selection
#   * structure-aware drag modes: soft, base, pair, loop, stem, branch, selection
#   * Alt-drag always moves only the clicked nucleotide
#   * smooth kinetic panning, centered zoom, fit/jump conveniences
#   * no Java, Jmol, browser, network, or third-party Python GUI dependency
# ============================================================================

__DSSR_PY2D_GEL_EDITOR_VERSION__ = "v4.0-gel-structure-editor"

import time as _dssr2d_time
from collections import deque as _dssr2d_deque


if QtWidgets is not None and QtCore is not None and QtGui is not None:
    _DSSR_PY2D_VIEW_BEFORE_V4 = Dssr2DGraphicsView
    _DSSR_PY2D_NODE_BEFORE_V4 = Dssr2DNodeItem
    _DSSR_PY2D_DIALOG_BEFORE_V4 = DssrPython2DDialog

    # ---------------------------------------------------------------------
    # Softer edge rendering: a faint halo under a crisp rounded line.
    # This is patched in place so the old edge constructor is not rebound.
    # ---------------------------------------------------------------------

    def _dssr2d_v4_edge_paint(self, painter, option, widget=None):
        try:
            painter.save()
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            kind = str(getattr(self, "kind", "backbone"))
            layer = int(getattr(self, "layer", 0))

            if kind == "backbone":
                color = QtGui.QColor(72, 88, 98, 170)
                width = 1.65
                style = QtCore.Qt.SolidLine
                halo_alpha = 22
            elif kind == "tertiary":
                color = QtGui.QColor(145, 88, 178, 145)
                width = 1.35
                style = QtCore.Qt.DashLine
                halo_alpha = 18
            elif layer > 0:
                palette = [
                    QtGui.QColor(226, 75, 84, 215),
                    QtGui.QColor(238, 139, 47, 215),
                    QtGui.QColor(174, 82, 199, 215),
                    QtGui.QColor(42, 158, 164, 215),
                ]
                color = palette[(layer - 1) % len(palette)]
                width = 2.05
                style = QtCore.Qt.SolidLine
                halo_alpha = 28
            else:
                color = QtGui.QColor(69, 112, 216, 220)
                width = 2.05
                style = QtCore.Qt.SolidLine
                halo_alpha = 30

            halo_color = QtGui.QColor(color)
            halo_color.setAlpha(halo_alpha)
            halo = QtGui.QPen(halo_color)
            halo.setWidthF(width + 4.2)
            halo.setStyle(style)
            halo.setCapStyle(QtCore.Qt.RoundCap)
            halo.setJoinStyle(QtCore.Qt.RoundJoin)
            painter.setPen(halo)
            painter.setBrush(QtCore.Qt.NoBrush)
            painter.drawPath(self.path())

            main = QtGui.QPen(color)
            main.setWidthF(width)
            main.setStyle(style)
            main.setCapStyle(QtCore.Qt.RoundCap)
            main.setJoinStyle(QtCore.Qt.RoundJoin)
            painter.setPen(main)
            painter.drawPath(self.path())
            painter.restore()
        except Exception:
            try:
                painter.restore()
            except Exception:
                pass
            try:
                QtWidgets.QGraphicsPathItem.paint(self, painter, option, widget)
            except Exception:
                pass

    try:
        Dssr2DEdgeItem.paint = _dssr2d_v4_edge_paint
        Dssr2DEdgeItem._gel_editor_v4_installed = True
    except Exception:
        pass


    class Dssr2DGraphicsViewV4(_DSSR_PY2D_VIEW_BEFORE_V4):
        """Easy canvas navigation: background drag pans, modifiers box-select."""

        def __init__(self, scene, parent=None):
            super(Dssr2DGraphicsViewV4, self).__init__(scene, parent)
            self.editor = parent
            self._v4_canvas_panning = False
            self._v4_pan_last = None
            self._v4_pan_samples = []
            self._v4_rubber = False
            self._v4_velocity = QtCore.QPointF(0.0, 0.0)
            self._v4_inertia_timer = QtCore.QTimer(self)
            self._v4_inertia_timer.setInterval(16)
            self._v4_inertia_timer.timeout.connect(self._v4_inertia_tick)
            try:
                self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
                self.setInteractive(True)
                self.setRubberBandSelectionMode(QtCore.Qt.IntersectsItemShape)
                self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
                self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorViewCenter)
                self.setViewportUpdateMode(
                    QtWidgets.QGraphicsView.BoundingRectViewportUpdate
                )
                self.setFocusPolicy(QtCore.Qt.StrongFocus)
            except Exception:
                pass

        @staticmethod
        def _v4_is_node_item(item):
            current = item
            for _ in range(6):
                if current is None:
                    return False
                if isinstance(current, Dssr2DNodeItemV4):
                    return True
                try:
                    current = current.parentItem()
                except Exception:
                    current = None
            return False

        def _v4_stop_inertia(self):
            try:
                self._v4_inertia_timer.stop()
            except Exception:
                pass
            self._v4_velocity = QtCore.QPointF(0.0, 0.0)

        def _v4_begin_pan(self, event):
            self._v4_stop_inertia()
            self._v4_canvas_panning = True
            self._v4_pan_last = event.pos()
            self._v4_pan_samples = []
            try:
                self.setCursor(QtCore.Qt.ClosedHandCursor)
            except Exception:
                pass
            event.accept()

        def mousePressEvent(self, event):
            try:
                item = self.itemAt(event.pos())
                button = event.button()
                modifiers = event.modifiers()
                on_node = self._v4_is_node_item(item)
            except Exception:
                item = None
                button = None
                modifiers = 0
                on_node = False

            if on_node:
                self._v4_stop_inertia()
                super(Dssr2DGraphicsViewV4, self).mousePressEvent(event)
                return

            select_modifier = bool(
                modifiers
                & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
            )
            if button == QtCore.Qt.LeftButton and select_modifier:
                self._v4_stop_inertia()
                self._v4_rubber = True
                try:
                    self.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
                except Exception:
                    pass
                super(Dssr2DGraphicsViewV4, self).mousePressEvent(event)
                return

            if button in (QtCore.Qt.LeftButton, QtCore.Qt.RightButton, QtCore.Qt.MiddleButton):
                self._v4_begin_pan(event)
                return

            super(Dssr2DGraphicsViewV4, self).mousePressEvent(event)

        def mouseMoveEvent(self, event):
            if self._v4_canvas_panning and self._v4_pan_last is not None:
                delta = event.pos() - self._v4_pan_last
                self._v4_pan_last = event.pos()
                now = _dssr2d_time.monotonic()
                self._v4_pan_samples.append((now, float(delta.x()), float(delta.y())))
                self._v4_pan_samples = self._v4_pan_samples[-8:]
                try:
                    self.horizontalScrollBar().setValue(
                        self.horizontalScrollBar().value() - delta.x()
                    )
                    self.verticalScrollBar().setValue(
                        self.verticalScrollBar().value() - delta.y()
                    )
                except Exception:
                    pass
                event.accept()
                return
            super(Dssr2DGraphicsViewV4, self).mouseMoveEvent(event)

        def mouseReleaseEvent(self, event):
            if self._v4_rubber:
                self._v4_rubber = False
                super(Dssr2DGraphicsViewV4, self).mouseReleaseEvent(event)
                try:
                    self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
                except Exception:
                    pass
                return

            if self._v4_canvas_panning:
                self._v4_canvas_panning = False
                self._v4_pan_last = None
                try:
                    self.unsetCursor()
                except Exception:
                    pass
                self._v4_start_inertia()
                event.accept()
                return
            super(Dssr2DGraphicsViewV4, self).mouseReleaseEvent(event)

        def _v4_start_inertia(self):
            samples = list(self._v4_pan_samples)
            self._v4_pan_samples = []
            if len(samples) < 2:
                return
            cutoff = samples[-1][0] - 0.12
            recent = [sample for sample in samples if sample[0] >= cutoff]
            if not recent:
                return
            dx = sum(sample[1] for sample in recent)
            dy = sum(sample[2] for sample in recent)
            elapsed = max(0.016, recent[-1][0] - recent[0][0])
            vx = max(-42.0, min(42.0, dx / elapsed * 0.016))
            vy = max(-42.0, min(42.0, dy / elapsed * 0.016))
            if math.hypot(vx, vy) < 1.5:
                return
            self._v4_velocity = QtCore.QPointF(vx, vy)
            self._v4_inertia_timer.start()

        def _v4_inertia_tick(self):
            vx = float(self._v4_velocity.x()) * 0.89
            vy = float(self._v4_velocity.y()) * 0.89
            self._v4_velocity = QtCore.QPointF(vx, vy)
            try:
                self.horizontalScrollBar().setValue(
                    self.horizontalScrollBar().value() - int(round(vx))
                )
                self.verticalScrollBar().setValue(
                    self.verticalScrollBar().value() - int(round(vy))
                )
            except Exception:
                self._v4_stop_inertia()
                return
            if math.hypot(vx, vy) < 0.35:
                self._v4_stop_inertia()

        def wheelEvent(self, event):
            try:
                delta = event.angleDelta().y()
            except Exception:
                try:
                    delta = event.delta()
                except Exception:
                    delta = 0
            if not delta:
                return
            current = abs(float(self.transform().m11()))
            factor = math.pow(1.00105, float(delta))
            target = current * factor
            if target < 0.08:
                factor = 0.08 / max(1.0e-9, current)
            elif target > 14.0:
                factor = 14.0 / max(1.0e-9, current)
            try:
                self.scale(factor, factor)
                event.accept()
            except Exception:
                super(Dssr2DGraphicsViewV4, self).wheelEvent(event)

        def keyPressEvent(self, event):
            key = event.key()
            if key == QtCore.Qt.Key_F:
                try:
                    self.editor.fit_scene()
                    event.accept()
                    return
                except Exception:
                    pass
            if key == QtCore.Qt.Key_C:
                try:
                    self.editor.fit_selected()
                    event.accept()
                    return
                except Exception:
                    pass
            mode_keys = {
                QtCore.Qt.Key_1: "soft",
                QtCore.Qt.Key_2: "base",
                QtCore.Qt.Key_3: "pair",
                QtCore.Qt.Key_4: "loop",
                QtCore.Qt.Key_5: "stem",
                QtCore.Qt.Key_6: "branch",
            }
            if key in mode_keys:
                try:
                    self.editor.set_drag_mode(mode_keys[key])
                    event.accept()
                    return
                except Exception:
                    pass
            super(Dssr2DGraphicsViewV4, self).keyPressEvent(event)

        def mouseDoubleClickEvent(self, event):
            try:
                item = self.itemAt(event.pos())
            except Exception:
                item = None
            if self._v4_is_node_item(item):
                current = item
                while current is not None and not isinstance(current, Dssr2DNodeItemV4):
                    try:
                        current = current.parentItem()
                    except Exception:
                        current = None
                if current is not None:
                    try:
                        self.editor.select_indices(
                            self.editor._v4_stem_group(current.nt_index),
                            replace=True,
                        )
                        self.editor.fit_selected()
                        event.accept()
                        return
                    except Exception:
                        pass
            try:
                self.editor.fit_scene()
                event.accept()
                return
            except Exception:
                pass
            super(Dssr2DGraphicsViewV4, self).mouseDoubleClickEvent(event)


    class Dssr2DNodeItemV4(_DSSR_PY2D_NODE_BEFORE_V4):
        """Glossy nucleotide with spring feedback and structure-aware dragging."""

        RADIUS = 14.0

        def __init__(self, viewer, nt, x, y):
            super(Dssr2DNodeItemV4, self).__init__(viewer, nt, x, y)
            self._v4_hover = False
            self._v4_pressed = False
            self._v4_visual_scale = 1.0
            self._v4_scale_target = 1.0
            self._v4_scale_velocity = 0.0
            self._v4_drag_weights = {}
            self._v4_last_drag_delta = QtCore.QPointF(0.0, 0.0)
            self._v4_last_move_pos = None
            self._v4_last_move_time = None
            self._v4_drag_speed = 0.0
            try:
                self.setAcceptHoverEvents(True)
                self.setCursor(QtCore.Qt.OpenHandCursor)
                self.setCacheMode(QtWidgets.QGraphicsItem.DeviceCoordinateCache)
            except Exception:
                pass

        def boundingRect(self):
            r = float(self.RADIUS)
            return QtCore.QRectF(-r - 7.0, -r - 7.0, 2.0 * r + 14.0, 2.0 * r + 14.0)

        def shape(self):
            path = QtGui.QPainterPath()
            r = float(self.RADIUS) + 3.5
            path.addEllipse(QtCore.QRectF(-r, -r, 2.0 * r, 2.0 * r))
            return path

        @staticmethod
        def _v4_alpha_color(color, alpha):
            result = QtGui.QColor(color)
            result.setAlpha(int(alpha))
            return result

        def paint(self, painter, option, widget=None):
            r = float(self.RADIUS)
            try:
                painter.save()
                painter.setRenderHint(QtGui.QPainter.Antialiasing, True)

                selected = bool(self.isSelected())
                hover = bool(self._v4_hover)
                pressed = bool(self._v4_pressed)
                base = QtGui.QColor(self._base_fill())

                # Soft shadow, drawn manually to remain fast while nodes move.
                shadow_alpha = 58 if (hover or selected) else 38
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(QtGui.QColor(20, 32, 46, shadow_alpha)))
                painter.drawEllipse(
                    QtCore.QRectF(-r + 2.6, -r + 4.0, 2.0 * r, 2.0 * r)
                )

                # Selection/hover aura.
                if selected or hover:
                    aura = QtGui.QColor(61, 164, 245, 82 if selected else 42)
                    painter.setPen(QtCore.Qt.NoPen)
                    painter.setBrush(QtGui.QBrush(aura))
                    aura_r = r + (4.6 if selected else 3.2)
                    painter.drawEllipse(
                        QtCore.QRectF(-aura_r, -aura_r, 2.0 * aura_r, 2.0 * aura_r)
                    )

                # Radial gel gradient: a bright upper-left core and deeper rim.
                gradient = QtGui.QRadialGradient(
                    QtCore.QPointF(-r * 0.38, -r * 0.46), r * 1.55
                )
                core = base.lighter(148)
                mid = base.lighter(112)
                rim = base.darker(112)
                core.setAlpha(252)
                mid.setAlpha(250)
                rim.setAlpha(248)
                gradient.setColorAt(0.00, QtGui.QColor(255, 255, 255, 248))
                gradient.setColorAt(0.20, core)
                gradient.setColorAt(0.66, mid)
                gradient.setColorAt(1.00, rim)

                border_color = QtGui.QColor(40, 55, 68, 205)
                border_width = 1.35
                if selected:
                    border_color = QtGui.QColor(39, 133, 231, 245)
                    border_width = 2.25
                elif hover:
                    border_color = QtGui.QColor(43, 151, 171, 235)
                    border_width = 1.85
                if pressed:
                    border_color = QtGui.QColor(26, 112, 208, 245)

                pen = QtGui.QPen(border_color)
                pen.setWidthF(border_width)
                pen.setCapStyle(QtCore.Qt.RoundCap)
                painter.setPen(pen)
                painter.setBrush(QtGui.QBrush(gradient))
                painter.drawEllipse(QtCore.QRectF(-r, -r, 2.0 * r, 2.0 * r))

                # Glassy highlight and a faint lower glow.
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 125)))
                painter.drawEllipse(
                    QtCore.QRectF(-r * 0.58, -r * 0.66, r * 0.92, r * 0.48)
                )
                painter.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 40)))
                painter.drawEllipse(
                    QtCore.QRectF(-r * 0.52, r * 0.18, r * 1.04, r * 0.48)
                )
                painter.restore()
            except Exception:
                try:
                    painter.restore()
                except Exception:
                    pass
                try:
                    QtWidgets.QGraphicsEllipseItem.paint(self, painter, option, widget)
                except Exception:
                    pass

        def _v4_set_target_scale(self, target, kick=0.0):
            self._v4_scale_target = float(target)
            self._v4_scale_velocity += float(kick)
            try:
                self.viewer._v4_ensure_animation()
            except Exception:
                pass

        def _v4_advance_visual(self):
            enabled = True
            try:
                enabled = bool(self.viewer.jelly_cb.isChecked())
            except Exception:
                pass
            target = self._v4_scale_target if enabled else 1.0
            stiffness = 0.24 if enabled else 0.42
            damping = 0.68 if enabled else 0.55
            self._v4_scale_velocity = (
                self._v4_scale_velocity
                + (target - self._v4_visual_scale) * stiffness
            ) * damping
            self._v4_visual_scale += self._v4_scale_velocity
            if (
                abs(target - self._v4_visual_scale) < 0.0006
                and abs(self._v4_scale_velocity) < 0.0006
            ):
                self._v4_visual_scale = target
                self._v4_scale_velocity = 0.0
            try:
                self.setScale(max(0.82, min(1.24, self._v4_visual_scale)))
                self.update()
            except Exception:
                pass
            return not (
                abs(target - self._v4_visual_scale) < 0.0007
                and abs(self._v4_scale_velocity) < 0.0007
            )

        def hoverEnterEvent(self, event):
            self._v4_hover = True
            self.setZValue(16.0)
            self._v4_set_target_scale(1.075, kick=0.018)
            self.update()
            try:
                super(Dssr2DNodeItemV4, self).hoverEnterEvent(event)
            except Exception:
                pass

        def hoverLeaveEvent(self, event):
            self._v4_hover = False
            self.setZValue(12.0 if self.isSelected() else 5.0)
            self._v4_set_target_scale(1.045 if self.isSelected() else 1.0)
            self.update()
            try:
                super(Dssr2DNodeItemV4, self).hoverLeaveEvent(event)
            except Exception:
                pass

        def itemChange(self, change, value):
            result = super(Dssr2DNodeItemV4, self).itemChange(change, value)
            try:
                if change == QtWidgets.QGraphicsItem.ItemSelectedHasChanged:
                    selected = bool(value)
                    self.setZValue(12.0 if selected else (16.0 if self._v4_hover else 5.0))
                    self._v4_set_target_scale(
                        1.045 if selected else (1.075 if self._v4_hover else 1.0),
                        kick=0.010 if selected else 0.0,
                    )
                    self.update()
            except Exception:
                pass
            return result

        def mousePressEvent(self, event):
            self._v4_pressed = True
            self._v4_set_target_scale(0.935, kick=-0.035)
            self._v4_last_move_pos = event.scenePos()
            self._v4_last_move_time = _dssr2d_time.monotonic()
            self._v4_drag_speed = 0.0
            super(Dssr2DNodeItemV4, self).mousePressEvent(event)
            try:
                if event.button() == QtCore.Qt.LeftButton and self._v3_dragging:
                    self.viewer._v4_prepare_node_drag(self, event.modifiers())
            except Exception:
                pass

        def mouseMoveEvent(self, event):
            if getattr(self, "_v3_pan_dragging", False):
                super(Dssr2DNodeItemV4, self).mouseMoveEvent(event)
                return
            if not getattr(self, "_v3_dragging", False) or self._v3_drag_origin is None:
                super(Dssr2DNodeItemV4, self).mouseMoveEvent(event)
                return

            delta = event.scenePos() - self._v3_drag_origin
            self._v4_last_drag_delta = QtCore.QPointF(delta)
            now = _dssr2d_time.monotonic()
            if self._v4_last_move_pos is not None and self._v4_last_move_time is not None:
                dt = max(0.001, now - self._v4_last_move_time)
                step = event.scenePos() - self._v4_last_move_pos
                instant = math.hypot(step.x(), step.y()) / dt
                self._v4_drag_speed = 0.72 * self._v4_drag_speed + 0.28 * instant
            self._v4_last_move_pos = event.scenePos()
            self._v4_last_move_time = now

            snap = False
            grid = 1.0
            try:
                snap = bool(self.viewer.snap_cb.isChecked())
                grid = max(1.0, float(self.viewer.grid_spin.value()))
            except Exception:
                pass

            jelly = True
            try:
                jelly = bool(self.viewer.jelly_cb.isChecked())
            except Exception:
                pass

            for index, start in self._v3_drag_starts.items():
                if index < 0 or index >= len(self.viewer.nodes):
                    continue
                weight = float(self._v4_drag_weights.get(index, 1.0))
                x = start.x() + delta.x() * weight
                y = start.y() + delta.y() * weight
                if snap:
                    x = round(x / grid) * grid
                    y = round(y / grid) * grid
                node = self.viewer.nodes[index]
                if jelly and weight < 0.999:
                    follow = 0.44 + 0.34 * weight
                    current = node.pos()
                    x = current.x() + (x - current.x()) * follow
                    y = current.y() + (y - current.y()) * follow
                node.setPos(x, y)
            event.accept()

        def mouseReleaseEvent(self, event):
            was_dragging = bool(getattr(self, "_v3_dragging", False))
            if was_dragging:
                # Finish followers at their intended weighted endpoint before the
                # v3 history snapshot is recorded.
                delta = QtCore.QPointF(self._v4_last_drag_delta)
                snap = False
                grid = 1.0
                try:
                    snap = bool(self.viewer.snap_cb.isChecked())
                    grid = max(1.0, float(self.viewer.grid_spin.value()))
                except Exception:
                    pass
                for index, start in self._v3_drag_starts.items():
                    if index < 0 or index >= len(self.viewer.nodes):
                        continue
                    weight = float(self._v4_drag_weights.get(index, 1.0))
                    x = start.x() + delta.x() * weight
                    y = start.y() + delta.y() * weight
                    if snap:
                        x = round(x / grid) * grid
                        y = round(y / grid) * grid
                    self.viewer.nodes[index].setPos(x, y)

            super(Dssr2DNodeItemV4, self).mouseReleaseEvent(event)
            self._v4_pressed = False
            target = 1.075 if self._v4_hover else (1.045 if self.isSelected() else 1.0)
            kick = min(0.11, max(0.02, self._v4_drag_speed * 0.00012)) if was_dragging else 0.025
            self._v4_set_target_scale(target, kick=kick)
            self._v4_drag_weights = {}
            self._v4_last_move_pos = None
            self._v4_last_move_time = None
            self.update()

        def contextMenuEvent(self, event):
            menu = QtWidgets.QMenu()
            header = menu.addAction(
                "%s  ·  nt %s"
                % (
                    str(self.nt.get("base", "N")),
                    str(self.nt.get("resi") or self.nt_index + 1),
                )
            )
            header.setEnabled(False)
            menu.addSeparator()

            actions = {}
            for key, text in (
                ("base", "Select base"),
                ("pair", "Select base pair"),
                ("loop", "Select loop / unpaired region"),
                ("stem", "Select stem"),
                ("branch", "Select whole branch"),
            ):
                actions[key] = menu.addAction(text)
            menu.addSeparator()
            center_action = menu.addAction("Center selection")
            reset_action = menu.addAction("Reset selection to automatic layout")
            undo_action = menu.addAction("Undo")

            try:
                chosen = menu.exec_(event.screenPos())
            except Exception:
                try:
                    chosen = menu.exec(event.screenPos())
                except Exception:
                    chosen = None

            if chosen == actions.get("base"):
                self.viewer.select_indices([self.nt_index], replace=True)
            elif chosen == actions.get("pair"):
                self.viewer.select_indices(
                    self.viewer._v4_pair_group(self.nt_index), replace=True
                )
            elif chosen == actions.get("loop"):
                self.viewer.select_indices(
                    self.viewer._v4_loop_group(self.nt_index), replace=True
                )
            elif chosen == actions.get("stem"):
                self.viewer.select_indices(
                    self.viewer._v4_stem_group(self.nt_index), replace=True
                )
            elif chosen == actions.get("branch"):
                self.viewer.select_indices(
                    self.viewer._v4_branch_group(self.nt_index), replace=True
                )
            elif chosen == center_action:
                self.viewer.fit_selected()
            elif chosen == reset_action:
                self.viewer.reset_selected_bases()
            elif chosen == undo_action:
                self.viewer.undo_layout()
            event.accept()


    # Public lookups used by the inherited scene builder now create v4 items.
    Dssr2DGraphicsView = Dssr2DGraphicsViewV4
    Dssr2DNodeItem = Dssr2DNodeItemV4


    class DssrPython2DDialogV4(_DSSR_PY2D_DIALOG_BEFORE_V4):
        """Gel-themed RNA editor with structure-aware drag modes."""

        def __init__(
            self,
            model,
            pymol_selection="all",
            algorithm="smart",
            number_every=10,
            show_tertiary=False,
            parent=None,
        ):
            self._v4_animating = False
            self._v4_adjacency = None
            super(DssrPython2DDialogV4, self).__init__(
                model=model,
                pymol_selection=pymol_selection,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=parent,
            )
            self._v4_install_gel_controls()
            self._v4_apply_visual_theme()
            self._v4_timer = QtCore.QTimer(self)
            self._v4_timer.setInterval(16)
            self._v4_timer.timeout.connect(self._v4_animation_tick)
            self._v4_ensure_animation()
            self._v3_update_editor_status("gel editor ready")

        def _add_edge(self, i, j, kind, layer=0, lw="", linear=False):
            edge = super(DssrPython2DDialogV4, self)._add_edge(
                i, j, kind, layer=layer, lw=lw, linear=linear
            )
            if edge is not None:
                try:
                    edge.setAcceptedMouseButtons(QtCore.Qt.NoButton)
                    edge.setAcceptHoverEvents(False)
                except Exception:
                    pass
            return edge

        def _v4_apply_visual_theme(self):
            try:
                self.setStyleSheet(
                    """
                    QDialog { background: #f4f8fb; }
                    QGroupBox {
                        border: 1px solid #c9d7e2;
                        border-radius: 10px;
                        margin-top: 8px;
                        padding-top: 8px;
                        background: rgba(255,255,255,205);
                        font-weight: 600;
                    }
                    QGroupBox::title {
                        subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px;
                        color: #344b5d;
                    }
                    QPushButton {
                        border: 1px solid #b9cbd8;
                        border-radius: 7px;
                        padding: 5px 10px;
                        background: #ffffff;
                        color: #243b4b;
                    }
                    QPushButton:hover { background: #e9f6ff; border-color: #78b8df; }
                    QPushButton:pressed { background: #d7efff; }
                    QPushButton:disabled { color: #9aa9b4; background: #eef2f5; }
                    QComboBox, QSpinBox, QDoubleSpinBox {
                        border: 1px solid #b9cbd8;
                        border-radius: 6px;
                        padding: 3px 6px;
                        background: white;
                    }
                    QCheckBox { spacing: 6px; }
                    """
                )
            except Exception:
                pass
            try:
                gradient = QtGui.QLinearGradient(0.0, 0.0, 0.0, 900.0)
                gradient.setColorAt(0.0, QtGui.QColor(255, 255, 255))
                gradient.setColorAt(0.55, QtGui.QColor(248, 252, 255))
                gradient.setColorAt(1.0, QtGui.QColor(235, 245, 251))
                self.view.setBackgroundBrush(QtGui.QBrush(gradient))
            except Exception:
                pass

        def _v4_install_gel_controls(self):
            group = QtWidgets.QGroupBox("gel interaction")
            row = QtWidgets.QHBoxLayout(group)

            row.addWidget(QtWidgets.QLabel("drag"))
            self.drag_mode_combo = QtWidgets.QComboBox()
            self.drag_mode_combo.addItem("soft neighborhood", "soft")
            self.drag_mode_combo.addItem("single base", "base")
            self.drag_mode_combo.addItem("base pair", "pair")
            self.drag_mode_combo.addItem("loop / unpaired", "loop")
            self.drag_mode_combo.addItem("stem", "stem")
            self.drag_mode_combo.addItem("whole branch", "branch")
            self.drag_mode_combo.addItem("current selection", "selection")
            self.drag_mode_combo.setCurrentIndex(0)
            self.drag_mode_combo.setToolTip(
                "1–6 switch modes. Alt-drag always moves only one base."
            )
            self.drag_mode_combo.currentIndexChanged.connect(
                lambda *_: self._v3_update_editor_status("drag mode changed")
            )
            row.addWidget(self.drag_mode_combo)

            self.jelly_cb = QtWidgets.QCheckBox("jelly motion")
            self.jelly_cb.setChecked(True)
            self.jelly_cb.setToolTip(
                "Spring hover/press feedback and elastic neighborhood following."
            )
            self.jelly_cb.toggled.connect(lambda *_: self._v4_ensure_animation())
            row.addWidget(self.jelly_cb)

            row.addWidget(QtWidgets.QLabel("follow"))
            self.follow_spin = QtWidgets.QDoubleSpinBox()
            self.follow_spin.setRange(0.10, 0.90)
            self.follow_spin.setSingleStep(0.05)
            self.follow_spin.setDecimals(2)
            self.follow_spin.setValue(0.55)
            self.follow_spin.setToolTip(
                "How much graph-neighboring nucleotides follow in soft mode."
            )
            row.addWidget(self.follow_spin)

            self.fit_selected_btn = QtWidgets.QPushButton("fit selected")
            self.fit_selected_btn.clicked.connect(self.fit_selected)
            row.addWidget(self.fit_selected_btn)

            row.addWidget(QtWidgets.QLabel("go to nt"))
            self.goto_spin = QtWidgets.QSpinBox()
            self.goto_spin.setRange(1, max(1, len(self.model.nts)))
            self.goto_spin.setValue(1)
            row.addWidget(self.goto_spin)
            self.goto_btn = QtWidgets.QPushButton("go")
            self.goto_btn.clicked.connect(self.goto_nucleotide)
            row.addWidget(self.goto_btn)

            hint = QtWidgets.QLabel(
                "drag base=edit · drag empty canvas=pan · Shift/Ctrl+drag canvas=box select · "
                "double-click base=stem · wheel=zoom · F=fit · C=fit selected"
            )
            hint.setWordWrap(True)
            hint.setStyleSheet("color:#557080; padding-left:6px;")
            row.addWidget(hint, 1)

            root = self.layout()
            if root is not None:
                try:
                    root.insertWidget(2, group)
                except Exception:
                    root.addWidget(group)
            self.gel_group = group

            # The old help label emphasizes middle-button panning. Replace its
            # text if present; no old widget is deleted.
            try:
                for label in self.editor_group.findChildren(QtWidgets.QLabel):
                    if "middle" in label.text().lower() or "space" in label.text().lower():
                        label.setText(
                            "drag bases directly · drag empty canvas to pan · "
                            "Shift/Ctrl+drag canvas for box selection · arrows nudge"
                        )
            except Exception:
                pass

        def drag_mode(self):
            try:
                data = self.drag_mode_combo.currentData()
                if data:
                    return str(data)
                return str(self.drag_mode_combo.currentText()).split()[0].lower()
            except Exception:
                return "soft"

        def set_drag_mode(self, mode):
            mode = str(mode or "soft").strip().lower()
            try:
                for index in range(self.drag_mode_combo.count()):
                    if str(self.drag_mode_combo.itemData(index)) == mode:
                        self.drag_mode_combo.setCurrentIndex(index)
                        return
            except Exception:
                pass

        def _v4_pair_table(self):
            return Dssr2DLayout._planar_pair_table(self.model)

        def _v4_pair_group(self, index):
            partner = self._v3_pair_partner(index)
            return [index] if partner < 0 else sorted({index, partner})

        def _v4_stem_group(self, index):
            return self._v3_stem_indices(index)

        def _v4_loop_group(self, index):
            n = len(self.model.nts)
            if index < 0 or index >= n:
                return []
            table = self._v4_pair_table()
            if table[index] >= 0:
                return self._v4_pair_group(index)

            left = index
            right = index
            while left > 0 and (left - 1) not in self.model.chain_breaks and table[left - 1] < 0:
                left -= 1
            while right + 1 < n and right not in self.model.chain_breaks and table[right + 1] < 0:
                right += 1
            result = set(range(left, right + 1))

            # Include a common closing pair when this is a hairpin/internal loop.
            if left > 0 and right + 1 < n:
                a = left - 1
                b = right + 1
                if table[a] == b:
                    result.add(a)
                    result.add(b)
            return sorted(result)

        def _v4_chain_segment(self, index):
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

        def _v4_branch_group(self, index):
            table = self._v4_pair_table()
            stems = Dssr2DLayout._stem_tree(table, self.model.chain_breaks)
            candidates = []
            for stem in stems:
                outer_i = int(stem.get("outer_i", -1))
                outer_j = int(stem.get("outer_j", -1))
                if outer_i <= index <= outer_j:
                    candidates.append(stem)
            if not candidates:
                return self._v4_chain_segment(index)
            stem = min(
                candidates,
                key=lambda item: int(item.get("outer_j", 0)) - int(item.get("outer_i", 0)),
            )
            return list(range(int(stem["outer_i"]), int(stem["outer_j"]) + 1))

        def _v4_build_adjacency(self):
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
            self._v4_adjacency = adjacency
            return adjacency

        def _v4_soft_weights(self, anchors):
            anchors = sorted({int(i) for i in anchors if 0 <= int(i) < len(self.model.nts)})
            if not anchors:
                return {}
            adjacency = self._v4_adjacency or self._v4_build_adjacency()
            strength = 0.55
            try:
                strength = float(self.follow_spin.value())
            except Exception:
                pass
            max_depth = 3
            distance = {index: 0 for index in anchors}
            queue = _dssr2d_deque(anchors)
            while queue:
                current = queue.popleft()
                depth = distance[current]
                if depth >= max_depth:
                    continue
                for neighbor in adjacency[current]:
                    if neighbor not in distance:
                        distance[neighbor] = depth + 1
                        queue.append(neighbor)
            weights = {}
            for index, depth in distance.items():
                weights[index] = 1.0 if depth == 0 else max(0.08, strength ** depth)
            return weights

        def _v4_prepare_node_drag(self, node, modifiers):
            index = int(node.nt_index)
            selected = sorted(
                item.nt_index for item in self.nodes if item.isSelected()
            )
            alt = bool(modifiers & QtCore.Qt.AltModifier)
            additive = bool(
                modifiers & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
            )
            mode = "base" if alt else self.drag_mode()

            if additive:
                group = selected or [index]
                weights = {item: 1.0 for item in group}
            elif mode == "pair":
                group = self._v4_pair_group(index)
                self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}
            elif mode == "loop":
                group = self._v4_loop_group(index)
                self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}
            elif mode == "stem":
                group = self._v4_stem_group(index)
                self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}
            elif mode == "branch":
                group = self._v4_branch_group(index)
                self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}
            elif mode == "selection":
                group = selected or [index]
                if not selected:
                    self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}
            elif mode == "soft":
                anchors = selected if (index in selected and len(selected) > 1) else [index]
                if anchors == [index]:
                    self.select_indices([index], replace=True)
                weights = self._v4_soft_weights(anchors)
                group = sorted(weights)
            else:  # base
                if index in selected and len(selected) > 1:
                    group = selected
                else:
                    group = [index]
                    self.select_indices(group, replace=True)
                weights = {item: 1.0 for item in group}

            node._v3_drag_starts = {
                item: QtCore.QPointF(self.nodes[item].pos())
                for item in group
                if 0 <= item < len(self.nodes)
            }
            node._v4_drag_weights = {
                item: float(weights.get(item, 1.0)) for item in node._v3_drag_starts
            }
            self._v3_update_editor_status("dragging %s" % mode)

        def fit_selected(self):
            selected = [node for node in self.nodes if node.isSelected()]
            if not selected:
                self.fit_scene()
                return
            rect = None
            for node in selected:
                current = node.sceneBoundingRect()
                rect = QtCore.QRectF(current) if rect is None else rect.united(current)
            if rect is None:
                return
            rect = rect.adjusted(-80.0, -80.0, 80.0, 80.0)
            try:
                self.view.fitInView(rect, QtCore.Qt.KeepAspectRatio)
            except Exception:
                pass

        def goto_nucleotide(self):
            index = int(self.goto_spin.value()) - 1
            if index < 0 or index >= len(self.nodes):
                return
            self.select_indices([index], replace=True)
            try:
                self.view.centerOn(self.nodes[index])
                current = abs(float(self.view.transform().m11()))
                if current < 1.4:
                    factor = 1.4 / max(1.0e-9, current)
                    self.view.scale(factor, factor)
            except Exception:
                pass
            self._v3_update_editor_status("centered nt %d" % (index + 1))

        def _v4_ensure_animation(self):
            self._v4_animating = True
            try:
                if hasattr(self, "_v4_timer") and not self._v4_timer.isActive():
                    self._v4_timer.start()
            except Exception:
                pass

        def _v4_animation_tick(self):
            active = False
            for node in list(getattr(self, "nodes", [])):
                try:
                    active = node._v4_advance_visual() or active
                except Exception:
                    pass
            if not active:
                self._v4_animating = False
                try:
                    self._v4_timer.stop()
                except Exception:
                    pass

        def redraw(self, *_args):
            super(DssrPython2DDialogV4, self).redraw(*_args)
            try:
                self._v4_adjacency = None
                self._v4_apply_visual_theme()
                self._v4_ensure_animation()
            except Exception:
                pass

        def _v3_update_editor_status(self, action=""):
            selected = sum(1 for node in getattr(self, "nodes", []) if node.isSelected())
            variant = str(
                getattr(self.model, "_dssr2d_layout_variant", self.algorithm)
            )
            text = (
                self.model.summary()
                + " | layout=%s (%s)" % (self.algorithm, variant)
                + " | selected=%d" % selected
                + " | drag=%s" % self.drag_mode()
                + " | background drag=pan; Shift/Ctrl+background drag=box select; Alt=single base"
            )
            if action:
                text += " | %s" % action
            try:
                self.status_label.setText(text)
            except Exception:
                pass


    DssrPython2DDialog = DssrPython2DDialogV4
    DssrPython2DDialog._gel_editor_v4_installed = True


try:
    print(
        "Loaded DSSR RNA 2D gel editor %s: glossy bases, elastic drag, "
        "easy canvas pan, structure-aware modes"
        % __DSSR_PY2D_GEL_EDITOR_VERSION__
    )
except Exception:
    pass

# ============================================================================
# DSSR pure-Python RNA 2D editor v5.0 — Standard scientific layout
#
# APPEND-ONLY PATCH: every byte above is retained.  This patch deliberately
# removes the *runtime effect* of the v4 gel theme without deleting its source.
# It provides one coherent publication-style visual language and a bundled
# pure-Python NAView layout, while preserving free nucleotide movement,
# multi-selection, undo/redo, save/load, PNG/SVG export, and PyMOL selection.
#
# The NAView implementation below is a Python adaptation of the Apache-2.0
# implementation in ViennaRNA/fornac (commit
# 36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1), itself based on the NAView
# algorithm.  The full Apache License 2.0 text is embedded immediately below.
# ============================================================================

__DSSR_PY2D_STANDARD_EDITOR_VERSION__ = "v5.0.0-standard-naview-editor"
__DSSR_FORNAC_NOTICE__ = (
    "Python adaptation of ViennaRNA/fornac src/naview/naview.js; "
    "fornac authors Peter Kerpedjiev, Stefan Hammer, and Ronny Lorenz; "
    "licensed under Apache-2.0."
)


__DSSR_FORNAC_APACHE_LICENSE__ = r"""

                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""


# ---------------------------------------------------------------------------
# Bundled NAView geometry engine (1-based internal indexing, planar scaffold)
# ---------------------------------------------------------------------------

class _DSSRNaviewRegion:
    __slots__ = ("start1", "end1", "start2", "end2")

    def __init__(self):
        self.start1 = 0
        self.end1 = 0
        self.start2 = 0
        self.end2 = 0


class _DSSRNaviewBase:
    __slots__ = ("mate", "x", "y", "extracted", "region")

    def __init__(self):
        self.mate = 0
        self.x = 9999.0
        self.y = 9999.0
        self.extracted = False
        self.region = None


class _DSSRNaviewConnection:
    __slots__ = (
        "loop",
        "region",
        "start",
        "end",
        "xrad",
        "yrad",
        "angle",
        "extruded",
        "broken",
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
        self.broken = False


class _DSSRNaviewLoop:
    __slots__ = (
        "connections",
        "number",
        "depth",
        "mark",
        "x",
        "y",
        "radius",
    )

    def __init__(self):
        self.connections = []
        self.number = 0
        self.depth = 0
        self.mark = False
        self.x = 0.0
        self.y = 0.0
        self.radius = 0.0

    @property
    def nconnection(self):
        return len(self.connections)


class _DSSRNaview:
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

        self.bases = [_DSSRNaviewBase() for _ in range(self.nbase + 1)]
        self.regions = [_DSSRNaviewRegion() for _ in range(self.nbase + 1)]
        self._read_in_bases(pair_table)
        self._find_regions()

        self.loops = [_DSSRNaviewLoop() for _ in range(self.nbase + 1)]
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
        result.number = self.loop_count
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
                        raise RuntimeError(
                            "NAView loop construction invariant failed"
                        )

                    connection = _DSSRNaviewConnection()
                    connection.loop = child_loop
                    connection.region = region
                    if index == region.start1:
                        connection.start = region.start1
                        connection.end = region.end2
                    else:
                        connection.start = region.start2
                        connection.end = region.end1
                    result.connections.append(connection)

                    reverse = _DSSRNaviewConnection()
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
                numerator / denominator
                if denominator > 1.0e-14
                else minimum_radius
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
                        separation = (
                            following.angle - candidate.angle
                        ) % (2.0 * math.pi)
                        if separation > largest_angle:
                            largest_angle = separation
                            largest_index = index
                    connection_end = largest_index
                    connection_start = (largest_index + 1) % loop.nconnection
                    loop.connections[connection_end].broken = True
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
                                half_angle = math.asin(
                                    min(1.0, 1.0 / (2.0 * radius))
                                )
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
                                following_index = (
                                    current_index + 1
                                ) % loop.nconnection
                                connection = loop.connections[current_index]
                                following = loop.connections[following_index]
                                angle = (
                                    connection.angle + following.angle
                                ) / 2.0
                                if connection.angle > following.angle:
                                    angle -= math.pi
                                line_x = math.sin(angle)
                                line_y = -math.cos(angle)
                                separation = (
                                    following.angle - connection.angle
                                ) % (2.0 * math.pi)
                                multiplier = (
                                    2.0
                                    if connection.extruded
                                    and separation <= math.pi / 2.0
                                    else (
                                        1.5 if connection.extruded else 1.0
                                    )
                                )
                                self.bases[connection.end].x = (
                                    self.bases[following.start].x
                                    + multiplier * line_x
                                )
                                self.bases[connection.end].y = (
                                    self.bases[following.start].y
                                    + multiplier * line_y
                                )
                                self.bases[connection.start].x = (
                                    self.bases[connection.end].x
                                    + connection.yrad
                                )
                                self.bases[connection.start].y = (
                                    self.bases[connection.end].y
                                    - connection.xrad
                                )
                            else:
                                previous_index = (
                                    current_index - 1
                                ) % loop.nconnection
                                previous = loop.connections[previous_index]
                                connection = loop.connections[current_index]
                                angle = (
                                    previous.angle + connection.angle
                                ) / 2.0
                                if previous.angle > connection.angle:
                                    angle -= math.pi
                                line_x = -math.sin(angle)
                                line_y = math.cos(angle)
                                separation = (
                                    connection.angle - previous.angle
                                ) % (2.0 * math.pi)
                                multiplier = (
                                    2.0
                                    if previous.extruded
                                    and separation <= math.pi / 2.0
                                    else (1.5 if previous.extruded else 1.0)
                                )
                                self.bases[connection.start].x = (
                                    self.bases[previous.end].x
                                    + multiplier * line_x
                                )
                                self.bases[connection.start].y = (
                                    self.bases[previous.end].y
                                    + multiplier * line_y
                                )
                                self.bases[connection.end].x = (
                                    self.bases[connection.start].x
                                    - connection.yrad
                                )
                                self.bases[connection.end].y = (
                                    self.bases[connection.start].y
                                    + connection.xrad
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
                expected = (
                    following.angle - connection.angle
                ) % (2.0 * math.pi)
                if abs(sweep - expected) > math.pi:
                    if not connection.extruded and (
                        following.start - connection.end
                    ) != 1:
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
                                local_radius = radius_current + (
                                    radius_following - radius_current
                                ) * (angle - angle_current) / sweep
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
            sum_x += (
                self.bases[connection.start].x
                + self.bases[connection.end].x
            )
            sum_y += (
                self.bases[connection.start].y
                + self.bases[connection.end].y
            )
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
                self.bases[base_index].x = (
                    self.bases[start].x + delta_x * offset / float(length)
                )
                self.bases[base_index].y = (
                    self.bases[start].y + delta_y * offset / float(length)
                )
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


def _dssr2d_v5_pair_table(model, start=0, end=None):
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


def _dssr2d_v5_segments(model):
    total = len(model.nts)
    if total <= 0:
        return []
    result = []
    start = 0
    for break_after in sorted(model.chain_breaks):
        break_after = int(break_after)
        if start <= break_after < total:
            result.append((start, break_after))
            start = break_after + 1
    if start < total:
        result.append((start, total - 1))
    return result or [(0, total - 1)]


def _dssr2d_v5_center(points):
    if not points:
        return []
    center_x = sum(point[0] for point in points) / float(len(points))
    center_y = sum(point[1] for point in points) / float(len(points))
    return [(point[0] - center_x, point[1] - center_y) for point in points]


def _dssr2d_v5_rotate(points, angle):
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return [
        (
            point[0] * cosine - point[1] * sine,
            point[0] * sine + point[1] * cosine,
        )
        for point in points
    ]


def _dssr2d_v5_standardize_trna_orientation(model, points):
    """Orient tRNA like the DSSR/VARNA 1EHZ example.

    The acceptor stem points down, the anticodon arm points up, the D arm is
    left, and the T arm is right.  Detection is topological, never PDB-name
    based.  Coordinates remain NAView coordinates; only a rigid rotation and
    optional mirror are applied.
    """
    try:
        topology = _dssr2d_v3_detect_trna_topology(model)
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
        points = _dssr2d_v5_rotate(points, math.pi / 2.0 - current_angle)

    # Mirror only when the earlier (D) arm appears to the right of the later
    # (T) arm.  This makes the orientation reproducible across structures.
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

    return _dssr2d_v5_center(points), True


def _dssr2d_v5_naview_layout(model):
    """Return one coherent NAView-style scientific layout."""
    total = len(model.nts)
    if total <= 0:
        return []

    scaffold_pair_count = len(model.planar_secondary_pairs())
    if scaffold_pair_count <= 0:
        model._dssr2d_layout_variant = "unpaired circular"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "circular")

    # NAView uses a circular exterior loop.  For multi-chain or intermolecular
    # complexes the old component-aware layout is safer than inventing a false
    # backbone link between chains.
    if model.chain_count() != 1:
        model._dssr2d_layout_variant = "multi-chain radiate fallback"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "radiate")

    table = _dssr2d_v5_pair_table(model)
    try:
        points = _DSSRNaview().coordinates(table)
    except Exception as error:
        try:
            model.warnings.append("NAView fallback: %s" % str(error))
        except Exception:
            pass
        model._dssr2d_layout_variant = "radiate fallback"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "radiate")

    # At 1.65x, the smallest NAView junction separation for 1EHZ remains larger
    # than the classic node diameter while the overall drawing still fits well.
    scale = 1.65
    points = [(scale * point[0], scale * point[1]) for point in points]
    points, is_trna = _dssr2d_v5_standardize_trna_orientation(model, points)
    model._dssr2d_layout_variant = (
        "NAView standard tRNA cloverleaf" if is_trna else "NAView standard"
    )
    return _dssr2d_v5_center(points)


_DSSR_PY2D_COMPUTE_BEFORE_V5 = Dssr2DLayout.compute


def _dssr2d_v5_compute(model, algorithm):
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
        return _dssr2d_v5_naview_layout(model)
    if name in ("legacy radiate", "legacy", "radiate", "radial"):
        model._dssr2d_layout_variant = "legacy radiate"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "radiate")
    if name in ("circular", "circle"):
        model._dssr2d_layout_variant = "circular"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "circular")
    if name in ("linear", "line"):
        model._dssr2d_layout_variant = "linear"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "linear")
    if name in ("force", "spring", "graph"):
        model._dssr2d_layout_variant = "force"
        return _DSSR_PY2D_COMPUTE_BEFORE_V5(model, "force")
    return _dssr2d_v5_naview_layout(model)


Dssr2DLayout.naview = staticmethod(_dssr2d_v5_naview_layout)
Dssr2DLayout.standard = staticmethod(_dssr2d_v5_naview_layout)
Dssr2DLayout.compute = staticmethod(_dssr2d_v5_compute)
Dssr2DLayout._standard_editor_v5_installed = True


# ---------------------------------------------------------------------------
# One uniform VARNA-like visual language: no gel, no side-dependent styling
# ---------------------------------------------------------------------------

if QtWidgets is not None and QtCore is not None and QtGui is not None:
    _DSSR_PY2D_VIEW_BEFORE_V5 = Dssr2DGraphicsView
    _DSSR_PY2D_NODE_BEFORE_V5 = Dssr2DNodeItem
    _DSSR_PY2D_DIALOG_BEFORE_V5 = DssrPython2DDialog

    def _dssr2d_v5_edge_style(self):
        kind = str(getattr(self, "kind", "backbone"))
        layer = int(getattr(self, "layer", 0))
        if kind == "backbone":
            color = QtGui.QColor(105, 105, 105)
            width = 1.05
            style = QtCore.Qt.SolidLine
        elif kind == "tertiary":
            color = QtGui.QColor(145, 115, 165)
            width = 0.95
            style = QtCore.Qt.DashLine
        else:
            # Primary and pseudoknot base pairs deliberately share exactly the
            # same blue scientific style, matching the simple VARNA convention.
            color = QtGui.QColor(60, 85, 175)
            width = 1.25 if layer == 0 else 1.15
            style = QtCore.Qt.SolidLine
        pen = QtGui.QPen(color)
        pen.setWidthF(width)
        pen.setStyle(style)
        try:
            pen.setCapStyle(QtCore.Qt.RoundCap)
            pen.setJoinStyle(QtCore.Qt.RoundJoin)
        except Exception:
            pass
        self.setPen(pen)
        if getattr(self, "lw", ""):
            self.setToolTip("Base pair: %s" % self.lw)

    def _dssr2d_v5_edge_paint(self, painter, option, widget=None):
        try:
            painter.save()
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            painter.setPen(self.pen())
            painter.setBrush(QtCore.Qt.NoBrush)
            painter.drawPath(self.path())
            painter.restore()
        except Exception:
            try:
                painter.restore()
            except Exception:
                pass
            QtWidgets.QGraphicsPathItem.paint(self, painter, option, widget)

    def _dssr2d_v5_edge_geometry(self):
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
                QtCore.QPointF(
                    0.5 * (first.x() + second.x()), sign * height
                ),
                second,
            )
        elif self.kind == "tertiary":
            delta_x = second.x() - first.x()
            delta_y = second.y() - first.y()
            distance = math.hypot(delta_x, delta_y)
            if distance <= 1.0e-9:
                path.lineTo(second)
            else:
                normal_x = -delta_y / distance
                normal_y = delta_x / distance
                curvature = min(120.0, max(24.0, 0.14 * distance))
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
            # Backbone, canonical pairs, and pseudoknot pairs are straight,
            # exactly as in the simple DSSR/VARNA 1EHZ diagram.
            path.lineTo(second)
        self.setPath(path)

    try:
        Dssr2DEdgeItem._set_style = _dssr2d_v5_edge_style
        Dssr2DEdgeItem.paint = _dssr2d_v5_edge_paint
        Dssr2DEdgeItem.update_geometry = _dssr2d_v5_edge_geometry
        Dssr2DEdgeItem._standard_editor_v5_installed = True
    except Exception:
        pass


    class Dssr2DGraphicsViewV5(_DSSR_PY2D_VIEW_BEFORE_V5):
        """V4's easy pan/zoom, with simpler scientific-editor shortcuts."""

        def keyPressEvent(self, event):
            mode_keys = {
                QtCore.Qt.Key_1: "base",
                QtCore.Qt.Key_2: "selection",
                QtCore.Qt.Key_3: "pair",
                QtCore.Qt.Key_4: "loop",
                QtCore.Qt.Key_5: "stem",
                QtCore.Qt.Key_6: "branch",
            }
            if event.key() in mode_keys:
                try:
                    self.editor.set_drag_mode(mode_keys[event.key()])
                    event.accept()
                    return
                except Exception:
                    pass
            _DSSR_PY2D_VIEW_BEFORE_V5.keyPressEvent(self, event)


    class Dssr2DNodeItemV5(_DSSR_PY2D_NODE_BEFORE_V5):
        """Uniform flat nucleotide glyph with direct, unrestricted dragging."""

        RADIUS = 9.5

        def __init__(self, viewer, nt, x, y):
            _DSSR_PY2D_NODE_BEFORE_V5.__init__(self, viewer, nt, x, y)
            self._v4_hover = False
            self._v4_pressed = False
            try:
                self.setScale(1.0)
                self.setAcceptHoverEvents(True)
                self.setCursor(QtCore.Qt.OpenHandCursor)
                self.setCacheMode(QtWidgets.QGraphicsItem.DeviceCoordinateCache)
            except Exception:
                pass
            try:
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(7)
                font.setBold(False)
                self.base_text_item.setFont(font)
                rect = self.base_text_item.boundingRect()
                self.base_text_item.setPos(
                    -0.5 * rect.width(), -0.5 * rect.height()
                )
                self.base_text_item.setBrush(
                    QtGui.QBrush(QtGui.QColor(25, 25, 25))
                )
            except Exception:
                pass

        def boundingRect(self):
            radius = float(self.RADIUS)
            return QtCore.QRectF(
                -radius - 3.0,
                -radius - 3.0,
                2.0 * radius + 6.0,
                2.0 * radius + 6.0,
            )

        def shape(self):
            radius = float(self.RADIUS) + 2.0
            path = QtGui.QPainterPath()
            path.addEllipse(
                QtCore.QRectF(-radius, -radius, 2.0 * radius, 2.0 * radius)
            )
            return path

        def _base_fill(self):
            if not bool(getattr(self.viewer, "base_colors", False)):
                return QtGui.QColor(252, 252, 252)
            base = str(self.nt.get("base", "N")).upper()[:1]
            palette = {
                "A": QtGui.QColor(233, 246, 229),
                "C": QtGui.QColor(229, 239, 250),
                "G": QtGui.QColor(250, 243, 216),
                "U": QtGui.QColor(250, 230, 230),
                "T": QtGui.QColor(250, 230, 230),
                "I": QtGui.QColor(239, 233, 248),
            }
            return palette.get(base, QtGui.QColor(245, 245, 245))

        def paint(self, painter, option, widget=None):
            radius = float(self.RADIUS)
            try:
                painter.save()
                painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
                selected = bool(self.isSelected())
                hovered = bool(getattr(self, "_v4_hover", False))

                if selected:
                    fill = QtGui.QColor(225, 238, 253)
                    border = QtGui.QColor(35, 105, 185)
                    width = 2.0
                elif hovered:
                    fill = QtGui.QColor(self._base_fill()).lighter(102)
                    border = QtGui.QColor(70, 105, 145)
                    width = 1.35
                else:
                    fill = QtGui.QColor(self._base_fill())
                    border = QtGui.QColor(48, 48, 48)
                    width = 1.05

                pen = QtGui.QPen(border)
                pen.setWidthF(width)
                painter.setPen(pen)
                painter.setBrush(QtGui.QBrush(fill))
                painter.drawEllipse(
                    QtCore.QRectF(
                        -radius, -radius, 2.0 * radius, 2.0 * radius
                    )
                )
                painter.restore()
            except Exception:
                try:
                    painter.restore()
                except Exception:
                    pass
                QtWidgets.QGraphicsEllipseItem.paint(self, painter, option, widget)

        # Disable all gel scaling while retaining V4's robust drag machinery.
        def _v4_set_target_scale(self, target, kick=0.0):
            try:
                self.setScale(1.0)
            except Exception:
                pass

        def _v4_advance_visual(self):
            try:
                self.setScale(1.0)
            except Exception:
                pass
            return False

        def hoverEnterEvent(self, event):
            self._v4_hover = True
            self.setZValue(12.0 if self.isSelected() else 8.0)
            self.update()
            try:
                QtWidgets.QGraphicsEllipseItem.hoverEnterEvent(self, event)
            except Exception:
                pass

        def hoverLeaveEvent(self, event):
            self._v4_hover = False
            self.setZValue(12.0 if self.isSelected() else 5.0)
            self.update()
            try:
                QtWidgets.QGraphicsEllipseItem.hoverLeaveEvent(self, event)
            except Exception:
                pass

        def itemChange(self, change, value):
            result = _DSSR_PY2D_NODE_BEFORE_V5.itemChange(self, change, value)
            try:
                if change == QtWidgets.QGraphicsItem.ItemSelectedHasChanged:
                    self.setScale(1.0)
                    self.setZValue(12.0 if bool(value) else 5.0)
                    self.update()
            except Exception:
                pass
            return result


    Dssr2DGraphicsView = Dssr2DGraphicsViewV5
    Dssr2DNodeItem = Dssr2DNodeItemV5


    class DssrPython2DDialogV5(_DSSR_PY2D_DIALOG_BEFORE_V5):
        """Consistent publication-style viewer plus a straightforward editor."""

        def __init__(
            self,
            model,
            pymol_selection="all",
            algorithm="standard",
            number_every=10,
            show_tertiary=False,
            parent=None,
        ):
            _DSSR_PY2D_DIALOG_BEFORE_V5.__init__(
                self,
                model=model,
                pymol_selection=pymol_selection,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=parent,
            )

            self.setWindowTitle("RNA 2D editor — %s" % self.model.title)
            try:
                self.resize(1120, 820)
            except Exception:
                pass

            # Replace every automatic-layout choice with a small, coherent set.
            try:
                self.layout_combo.blockSignals(True)
                self.layout_combo.clear()
                self.layout_combo.addItems(
                    ["standard", "circular", "linear", "legacy radiate"]
                )
                selected = self.layout_combo.findText(
                    str(algorithm or "standard").strip().lower()
                )
                self.layout_combo.setCurrentIndex(selected if selected >= 0 else 0)
                self.layout_combo.blockSignals(False)
            except Exception:
                pass

            try:
                self.number_spin.blockSignals(True)
                self.number_spin.setValue(max(0, int(number_every)))
                self.number_spin.blockSignals(False)
                self.tertiary_cb.blockSignals(True)
                self.tertiary_cb.setChecked(bool(show_tertiary))
                self.tertiary_cb.blockSignals(False)
                self.base_colors_cb.blockSignals(True)
                self.base_colors_cb.setChecked(False)
                self.base_colors_cb.blockSignals(False)
            except Exception:
                pass

            self._v3_force_relayout = True
            self._v3_undo = []
            self._v3_redo = []
            self.redraw()
            self._v3_update_history_buttons()
            self._v3_update_editor_status("standard editor ready")
            try:
                QtCore.QTimer.singleShot(0, self.fit_scene)
            except Exception:
                pass

        def _v4_apply_visual_theme(self):
            try:
                self.setStyleSheet(
                    """
                    QDialog { background: #f2f2f2; }
                    QGroupBox {
                        border: 1px solid #c8c8c8;
                        border-radius: 3px;
                        margin-top: 7px;
                        padding-top: 7px;
                        background: #f7f7f7;
                    }
                    QGroupBox::title {
                        subcontrol-origin: margin;
                        left: 8px;
                        padding: 0 4px;
                        color: #333333;
                    }
                    QPushButton, QComboBox, QSpinBox, QDoubleSpinBox {
                        min-height: 22px;
                    }
                    """
                )
            except Exception:
                pass
            try:
                self.view.setBackgroundBrush(
                    QtGui.QBrush(QtGui.QColor(255, 255, 255))
                )
            except Exception:
                pass

        def _v4_install_gel_controls(self):
            # Method name is retained because V4 calls it during construction;
            # the resulting group is deliberately a plain editing group.
            group = QtWidgets.QGroupBox("editing")
            row = QtWidgets.QHBoxLayout(group)

            row.addWidget(QtWidgets.QLabel("drag"))
            self.drag_mode_combo = QtWidgets.QComboBox()
            self.drag_mode_combo.addItem("single base", "base")
            self.drag_mode_combo.addItem("current selection", "selection")
            self.drag_mode_combo.addItem("base pair", "pair")
            self.drag_mode_combo.addItem("loop / unpaired", "loop")
            self.drag_mode_combo.addItem("stem", "stem")
            self.drag_mode_combo.addItem("whole branch", "branch")
            self.drag_mode_combo.setCurrentIndex(0)
            self.drag_mode_combo.currentIndexChanged.connect(
                lambda *_args: self._v3_update_editor_status("drag mode changed")
            )
            row.addWidget(self.drag_mode_combo)

            self.fit_selected_btn = QtWidgets.QPushButton("fit selected")
            self.fit_selected_btn.clicked.connect(self.fit_selected)
            row.addWidget(self.fit_selected_btn)

            row.addWidget(QtWidgets.QLabel("go to nt"))
            self.goto_spin = QtWidgets.QSpinBox()
            self.goto_spin.setRange(1, max(1, len(self.model.nts)))
            self.goto_spin.setValue(1)
            row.addWidget(self.goto_spin)
            self.goto_btn = QtWidgets.QPushButton("go")
            self.goto_btn.clicked.connect(self.goto_nucleotide)
            row.addWidget(self.goto_btn)

            hint = QtWidgets.QLabel(
                "base drag=move · selected bases move together · empty drag=pan · "
                "Shift/Ctrl+empty drag=box select · wheel=zoom · Ctrl+Z/Y=undo/redo"
            )
            hint.setWordWrap(True)
            row.addWidget(hint, 1)

            # Compatibility attributes expected by inherited V4 methods.  They
            # are intentionally hidden and disabled.
            self.jelly_cb = QtWidgets.QCheckBox()
            self.jelly_cb.setChecked(False)
            self.jelly_cb.setVisible(False)
            self.follow_spin = QtWidgets.QDoubleSpinBox()
            self.follow_spin.setValue(0.0)
            self.follow_spin.setVisible(False)

            root = self.layout()
            if root is not None:
                try:
                    root.insertWidget(2, group)
                except Exception:
                    root.addWidget(group)
            self.gel_group = group

            try:
                for label in self.editor_group.findChildren(QtWidgets.QLabel):
                    if "middle" in label.text().lower() or "space" in label.text().lower():
                        label.setText(
                            "drag bases directly · Shift/Ctrl click for multi-select · "
                            "arrows nudge · save/load preserves manual coordinates"
                        )
            except Exception:
                pass

        def _v4_ensure_animation(self):
            return

        def _v4_animation_tick(self):
            try:
                if hasattr(self, "_v4_timer"):
                    self._v4_timer.stop()
            except Exception:
                pass

        def drag_mode(self):
            try:
                data = self.drag_mode_combo.currentData()
                return str(data or "base")
            except Exception:
                return "base"

        def set_drag_mode(self, mode):
            wanted = str(mode or "base").strip().lower()
            try:
                for index in range(self.drag_mode_combo.count()):
                    if str(self.drag_mode_combo.itemData(index)) == wanted:
                        self.drag_mode_combo.setCurrentIndex(index)
                        return
            except Exception:
                pass

        def redraw(self, *_args):
            _DSSR_PY2D_DIALOG_BEFORE_V5.redraw(self, *_args)
            try:
                self._v4_apply_visual_theme()
                for node in getattr(self, "nodes", []):
                    node.setScale(1.0)
            except Exception:
                pass

        def _add_number_labels(self):
            """Add sparse residue numbers without covering bases or earlier labels."""
            total = len(self.nodes)
            if total <= 0:
                return
            period = int(self.number_every)
            indices = {0, total - 1}
            if period > 0:
                indices.update(
                    index
                    for index in range(total)
                    if (index + 1) % period == 0
                )
            for break_after in self.model.chain_breaks:
                if 0 <= break_after < total:
                    indices.add(break_after)
                if 0 <= break_after + 1 < total:
                    indices.add(break_after + 1)

            center_x = sum(node.pos().x() for node in self.nodes) / float(total)
            center_y = sum(node.pos().y() for node in self.nodes) / float(total)
            node_radius = float(getattr(Dssr2DNodeItem, "RADIUS", 9.5)) + 2.5
            node_rects = [
                QtCore.QRectF(
                    node.pos().x() - node_radius,
                    node.pos().y() - node_radius,
                    2.0 * node_radius,
                    2.0 * node_radius,
                )
                for node in self.nodes
            ]
            used_label_rects = []

            def _intersection_area(first, second):
                try:
                    overlap = first.intersected(second)
                    if overlap.isEmpty():
                        return 0.0
                    return max(0.0, overlap.width()) * max(0.0, overlap.height())
                except Exception:
                    return 0.0

            for index in sorted(indices):
                node = self.nodes[index]
                nt = self.model.nts[index]
                value = str(nt.get("resi") or nt.get("number", index + 1))
                label = QtWidgets.QGraphicsSimpleTextItem(value, node)
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(7)
                label.setFont(font)
                label.setBrush(QtGui.QBrush(QtGui.QColor(35, 35, 35)))

                previous = None
                following = None
                if index > 0 and (index - 1) not in self.model.chain_breaks:
                    previous = self.nodes[index - 1].pos()
                if index + 1 < total and index not in self.model.chain_breaks:
                    following = self.nodes[index + 1].pos()

                if previous is not None and following is not None:
                    tangent_x = following.x() - previous.x()
                    tangent_y = following.y() - previous.y()
                elif following is not None:
                    tangent_x = following.x() - node.pos().x()
                    tangent_y = following.y() - node.pos().y()
                elif previous is not None:
                    tangent_x = node.pos().x() - previous.x()
                    tangent_y = node.pos().y() - previous.y()
                else:
                    tangent_x, tangent_y = 1.0, 0.0

                tangent_length = math.hypot(tangent_x, tangent_y)
                if tangent_length <= 1.0e-8:
                    tangent_x, tangent_y, tangent_length = 1.0, 0.0, 1.0
                tangent_x /= tangent_length
                tangent_y /= tangent_length
                normal_x, normal_y = -tangent_y, tangent_x

                outward_x = node.pos().x() - center_x
                outward_y = node.pos().y() - center_y
                outward_length = math.hypot(outward_x, outward_y)
                if outward_length <= 1.0e-8:
                    outward_x, outward_y, outward_length = normal_x, normal_y, 1.0
                outward_x /= outward_length
                outward_y /= outward_length
                if normal_x * outward_x + normal_y * outward_y < 0.0:
                    normal_x = -normal_x
                    normal_y = -normal_y

                directions = [
                    (normal_x, normal_y),
                    (outward_x, outward_y),
                    (-normal_x, -normal_y),
                    (tangent_x, tangent_y),
                    (-tangent_x, -tangent_y),
                ]
                distances = (17.0, 22.0, 28.0, 34.0)
                rect = label.boundingRect()
                best = None
                for direction_rank, (dir_x, dir_y) in enumerate(directions):
                    for distance_rank, distance in enumerate(distances):
                        local_x = dir_x * distance - 0.5 * rect.width()
                        local_y = dir_y * distance - 0.5 * rect.height()
                        scene_rect = QtCore.QRectF(
                            node.pos().x() + local_x,
                            node.pos().y() + local_y,
                            rect.width(),
                            rect.height(),
                        )
                        overlap = 0.0
                        for node_index, occupied in enumerate(node_rects):
                            if node_index == index:
                                continue
                            overlap += 8.0 * _intersection_area(scene_rect, occupied)
                        for occupied in used_label_rects:
                            overlap += 12.0 * _intersection_area(scene_rect, occupied)
                        # Prefer the outward local normal and a short leader-free offset
                        # whenever collision scores are equal.
                        score = overlap + 0.08 * distance_rank + 0.04 * direction_rank
                        candidate = (score, local_x, local_y, scene_rect)
                        if best is None or candidate[0] < best[0]:
                            best = candidate
                if best is None:
                    best = (0.0, 16.0, -16.0, QtCore.QRectF())
                label.setPos(best[1], best[2])
                used_label_rects.append(best[3])
                label.setZValue(8.0)
                try:
                    label.setAcceptedMouseButtons(QtCore.Qt.NoButton)
                except Exception:
                    pass

        def _v3_update_editor_status(self, action=""):
            selected = sum(
                1 for node in getattr(self, "nodes", []) if node.isSelected()
            )
            variant = str(
                getattr(self.model, "_dssr2d_layout_variant", self.algorithm)
            )
            text = (
                self.model.summary()
                + " | %s" % variant
                + " | selected=%d" % selected
                + " | drag=%s" % self.drag_mode()
                + " | left-drag base=move; left-drag empty=pan; Shift/Ctrl=multi-select"
            )
            if action:
                text += " | %s" % action
            try:
                self.status_label.setText(text)
            except Exception:
                pass


    DssrPython2DDialog = DssrPython2DDialogV5
    DssrPython2DDialog._standard_editor_v5_installed = True


# ---------------------------------------------------------------------------
# Main DSSR GUI and command defaults
# ---------------------------------------------------------------------------

_DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V5 = (
    DssrPython2DIntegration.install_gui_controls
)
_DSSR_PY2D_OPEN_DIALOG_BEFORE_V5 = DssrPython2DIntegration._open_dialog


def _dssr2d_v5_install_gui_controls(dialog):
    _DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V5(dialog)
    try:
        combo = dialog.py2d_layout_combo
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(["standard", "circular", "linear", "legacy radiate"])
        combo.setCurrentIndex(0)
        combo.blockSignals(False)
        dialog.py2d_number_spin.setValue(10)
        dialog.py2d_tertiary_cb.setChecked(False)
        dialog.py2d_group.setTitle(
            "RNA 2D (Python — standard scientific layout)"
        )
        dialog.py2d_open_btn.setToolTip(
            "Open a consistent NAView/VARNA-style editor. Drag any nucleotide "
            "directly; no Java, Jmol, ZIP, browser, or network is required."
        )
    except Exception:
        pass


def _dssr2d_v5_open_dialog(
    dssr_data,
    selection="all",
    title="RNA secondary structure",
    algorithm="standard",
    number_every=10,
    show_tertiary=0,
    parent=None,
):
    return _DSSR_PY2D_OPEN_DIALOG_BEFORE_V5(
        dssr_data=dssr_data,
        selection=selection,
        title=title,
        algorithm=algorithm,
        number_every=number_every,
        show_tertiary=show_tertiary,
        parent=parent,
    )


DssrPython2DIntegration.install_gui_controls = staticmethod(
    _dssr2d_v5_install_gui_controls
)
DssrPython2DIntegration._open_dialog = staticmethod(_dssr2d_v5_open_dialog)
DssrPython2DIntegration._standard_editor_v5_installed = True


_DSSR_PY2D_COMMAND_BEFORE_V5 = dssr_2d


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
    """Open the standard pure-Python RNA 2D editor.

    Layout choices are ``standard`` (bundled NAView), ``circular``, ``linear``,
    and ``legacy radiate``.  The default visual style follows the simple
    DSSR/VARNA convention: white circular residues, gray backbone, blue base
    pairs, numbering every ten residues, and hidden optional tertiary contacts.
    """
    return _DSSR_PY2D_COMMAND_BEFORE_V5(
        selection=selection,
        state=state,
        exe=exe,
        layout=layout,
        number_every=number_every,
        show_tertiary=show_tertiary,
        title=title,
        quiet=quiet,
    )


DssrFunctions.dssr_2d = staticmethod(dssr_2d)
cmd.extend("dssr_2d", dssr_2d)

try:
    print(
        "Loaded DSSR RNA 2D standard editor %s: bundled NAView, uniform "
        "VARNA-like styling, direct base dragging"
        % __DSSR_PY2D_STANDARD_EDITOR_VERSION__
    )
except Exception:
    pass

# ============================================================================
# DSSR pure-Python RNA 2D editor v6.0 — Premium gel + free brush selection
#
# APPEND-ONLY PATCH: retain the v5 scientific editor and add a reversible,
# coloured gel presentation plus a paint-like selection tool.  The selected
# 2D nucleotides continue to drive PyMOL's ``sele`` selection and are also
# rendered as a separate, non-destructive 3D highlight object.
# ============================================================================

__DSSR_PY2D_PREMIUM_EDITOR_VERSION__ = "v6.1.2-white-canvas"


if QtWidgets is not None and QtCore is not None and QtGui is not None:
    _DSSR_PY2D_VIEW_BEFORE_V6 = Dssr2DGraphicsView
    _DSSR_PY2D_NODE_BEFORE_V6 = Dssr2DNodeItem
    _DSSR_PY2D_DIALOG_BEFORE_V6 = DssrPython2DDialog
    _DSSR_PY2D_EDGE_PAINT_BEFORE_V6 = Dssr2DEdgeItem.paint

    def _dssr2d_v6_edge_paint(self, painter, option, widget=None):
        """Paint luminous gel edges, or delegate to v5 in scientific mode."""
        gel = False
        try:
            gel = bool(self.node_a.viewer.gel_style_enabled())
        except Exception:
            pass
        if not gel:
            return _DSSR_PY2D_EDGE_PAINT_BEFORE_V6(
                self, painter, option, widget
            )

        kind = str(getattr(self, "kind", "backbone"))
        layer = int(getattr(self, "layer", 0))
        if kind == "backbone":
            color = QtGui.QColor(112, 154, 186, 225)
            width = 1.65
            style = QtCore.Qt.SolidLine
        elif kind == "tertiary":
            color = QtGui.QColor(235, 111, 255, 220)
            width = 1.55
            style = QtCore.Qt.DashLine
        elif layer > 0:
            palette = (
                QtGui.QColor(183, 123, 255, 240),
                QtGui.QColor(255, 112, 176, 240),
                QtGui.QColor(255, 190, 82, 240),
            )
            color = palette[(layer - 1) % len(palette)]
            width = 2.05
            style = QtCore.Qt.SolidLine
        else:
            color = QtGui.QColor(78, 200, 255, 245)
            width = 2.15
            style = QtCore.Qt.SolidLine

        saved = False
        try:
            painter.save()
            saved = True
            painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
            glow_color = QtGui.QColor(color)
            glow_color.setAlpha(48)
            glow = QtGui.QPen(glow_color)
            glow.setWidthF(width + 4.2)
            glow.setStyle(style)
            glow.setCapStyle(QtCore.Qt.RoundCap)
            glow.setJoinStyle(QtCore.Qt.RoundJoin)
            painter.setPen(glow)
            painter.setBrush(QtCore.Qt.NoBrush)
            painter.drawPath(self.path())

            main = QtGui.QPen(color)
            main.setWidthF(width)
            main.setStyle(style)
            main.setCapStyle(QtCore.Qt.RoundCap)
            main.setJoinStyle(QtCore.Qt.RoundJoin)
            painter.setPen(main)
            painter.drawPath(self.path())
            painter.restore()
        except Exception:
            if saved:
                try:
                    painter.restore()
                except Exception:
                    pass
            return _DSSR_PY2D_EDGE_PAINT_BEFORE_V6(
                self, painter, option, widget
            )

    Dssr2DEdgeItem.paint = _dssr2d_v6_edge_paint
    Dssr2DEdgeItem._premium_editor_v6_installed = True


    class Dssr2DGraphicsViewV6(_DSSR_PY2D_VIEW_BEFORE_V6):
        """Canvas navigation plus a freehand, radius-based nucleotide brush."""

        def __init__(self, scene, parent=None):
            _DSSR_PY2D_VIEW_BEFORE_V6.__init__(self, scene, parent)
            self.editor = parent
            self._v6_brushing = False
            self._v6_brush_erase = False
            self._v6_brush_last = None
            self._v6_brush_path = None
            self._v6_brush_item = None
            self._v6_brush_changed = False

        def _v6_tool(self):
            try:
                return self.editor.interaction_tool()
            except Exception:
                return "edit"

        def _v6_scene_radius(self):
            radius_px = 34.0
            try:
                radius_px = float(self.editor.brush_radius_spin.value())
            except Exception:
                pass
            try:
                scale = abs(float(self.transform().m11()))
            except Exception:
                scale = 1.0
            return radius_px / max(0.08, scale)

        @staticmethod
        def _v6_distance_to_segment(point, first, second):
            px = float(point.x())
            py = float(point.y())
            ax = float(first.x())
            ay = float(first.y())
            bx = float(second.x())
            by = float(second.y())
            dx = bx - ax
            dy = by - ay
            length_sq = dx * dx + dy * dy
            if length_sq <= 1.0e-12:
                return math.hypot(px - ax, py - ay)
            ratio = ((px - ax) * dx + (py - ay) * dy) / length_sq
            ratio = max(0.0, min(1.0, ratio))
            return math.hypot(px - (ax + ratio * dx), py - (ay + ratio * dy))

        def _v6_begin_brush(self, event):
            self._v4_stop_inertia()
            self._v6_brushing = True
            self._v6_brush_erase = bool(event.modifiers() & QtCore.Qt.AltModifier)
            additive = bool(
                event.modifiers()
                & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier)
            )
            point = self.mapToScene(event.pos())
            self._v6_brush_last = QtCore.QPointF(point)
            self._v6_brush_changed = False

            if not additive and not self._v6_brush_erase:
                try:
                    self.editor._v3_rebuilding = True
                    self.scene().clearSelection()
                finally:
                    self.editor._v3_rebuilding = False

            path = QtGui.QPainterPath(point)
            self._v6_brush_path = path
            try:
                item = QtWidgets.QGraphicsPathItem()
                item.setPath(path)
                color = QtGui.QColor(255, 100, 176, 70)
                if not self._v6_brush_erase:
                    color = QtGui.QColor(58, 214, 255, 68)
                pen = QtGui.QPen(color)
                pen.setWidthF(2.0 * self._v6_scene_radius())
                pen.setCapStyle(QtCore.Qt.RoundCap)
                pen.setJoinStyle(QtCore.Qt.RoundJoin)
                item.setPen(pen)
                item.setBrush(QtCore.Qt.NoBrush)
                item.setZValue(3.0)
                item.setAcceptedMouseButtons(QtCore.Qt.NoButton)
                self.scene().addItem(item)
                self._v6_brush_item = item
            except Exception:
                self._v6_brush_item = None

            self._v6_apply_brush_segment(point, point)
            try:
                self.setCursor(QtCore.Qt.CrossCursor)
            except Exception:
                pass
            event.accept()

        def _v6_apply_brush_segment(self, first, second):
            radius = self._v6_scene_radius()
            changed = False
            try:
                self.editor._v3_rebuilding = True
                for node in list(getattr(self.editor, "nodes", [])):
                    node_radius = float(getattr(node, "RADIUS", 10.0)) * 0.65
                    distance = self._v6_distance_to_segment(
                        node.scenePos(), first, second
                    )
                    if distance <= radius + node_radius:
                        wanted = not self._v6_brush_erase
                        if bool(node.isSelected()) != wanted:
                            node.setSelected(wanted)
                            changed = True
            finally:
                self.editor._v3_rebuilding = False
            self._v6_brush_changed = self._v6_brush_changed or changed
            if changed:
                self.editor._v6_after_brush_selection(final=False)

        def _v6_finish_brush(self, event=None):
            if not self._v6_brushing:
                return
            self._v6_brushing = False
            self._v6_brush_last = None
            if self._v6_brush_item is not None:
                try:
                    self.scene().removeItem(self._v6_brush_item)
                except Exception:
                    pass
            self._v6_brush_item = None
            self._v6_brush_path = None
            self.editor._v6_after_brush_selection(final=True)
            try:
                self.editor._v6_update_view_cursor()
            except Exception:
                pass
            if event is not None:
                event.accept()

        def mousePressEvent(self, event):
            if (
                event.button() == QtCore.Qt.LeftButton
                and self._v6_tool() == "brush"
            ):
                self._v6_begin_brush(event)
                return
            _DSSR_PY2D_VIEW_BEFORE_V6.mousePressEvent(self, event)

        def mouseMoveEvent(self, event):
            if self._v6_brushing and self._v6_brush_last is not None:
                point = self.mapToScene(event.pos())
                first = QtCore.QPointF(self._v6_brush_last)
                self._v6_brush_last = QtCore.QPointF(point)
                if self._v6_brush_path is not None:
                    self._v6_brush_path.lineTo(point)
                    if self._v6_brush_item is not None:
                        self._v6_brush_item.setPath(self._v6_brush_path)
                self._v6_apply_brush_segment(first, point)
                event.accept()
                return
            _DSSR_PY2D_VIEW_BEFORE_V6.mouseMoveEvent(self, event)

        def mouseReleaseEvent(self, event):
            if self._v6_brushing and event.button() == QtCore.Qt.LeftButton:
                self._v6_finish_brush(event)
                return
            _DSSR_PY2D_VIEW_BEFORE_V6.mouseReleaseEvent(self, event)

        def keyPressEvent(self, event):
            if event.key() == QtCore.Qt.Key_B:
                self.editor.set_interaction_tool("brush")
                event.accept()
                return
            if event.key() == QtCore.Qt.Key_P:
                self.editor.set_interaction_tool("edit")
                event.accept()
                return
            _DSSR_PY2D_VIEW_BEFORE_V6.keyPressEvent(self, event)


    class Dssr2DNodeItemV6(_DSSR_PY2D_NODE_BEFORE_V6):
        """A coloured nucleotide that switches between gel and flat rendering."""

        RADIUS = 12.5

        def __init__(self, viewer, nt, x, y):
            _DSSR_PY2D_NODE_BEFORE_V6.__init__(self, viewer, nt, x, y)
            try:
                font = QtGui.QFont("Sans Serif")
                font.setPointSize(8)
                font.setBold(True)
                self.base_text_item.setFont(font)
                rect = self.base_text_item.boundingRect()
                self.base_text_item.setPos(-0.5 * rect.width(), -0.5 * rect.height())
            except Exception:
                pass

        def boundingRect(self):
            radius = float(self.RADIUS)
            return QtCore.QRectF(
                -radius - 7.0, -radius - 7.0,
                2.0 * radius + 14.0, 2.0 * radius + 14.0,
            )

        def shape(self):
            radius = float(self.RADIUS) + 3.0
            path = QtGui.QPainterPath()
            path.addEllipse(
                QtCore.QRectF(-radius, -radius, 2.0 * radius, 2.0 * radius)
            )
            return path

        def _base_fill(self):
            use_colors = bool(getattr(self.viewer, "base_colors", True))
            if not use_colors:
                if self.viewer.gel_style_enabled():
                    return QtGui.QColor(184, 207, 222)
                return QtGui.QColor(250, 252, 254)
            base = str(self.nt.get("base", "N")).upper()[:1]
            if self.viewer.gel_style_enabled():
                palette = {
                    "A": QtGui.QColor(255, 91, 134),
                    "C": QtGui.QColor(57, 169, 255),
                    "G": QtGui.QColor(255, 191, 66),
                    "U": QtGui.QColor(48, 211, 171),
                    "T": QtGui.QColor(164, 118, 255),
                    "I": QtGui.QColor(191, 116, 255),
                }
                return palette.get(base, QtGui.QColor(132, 151, 176))
            palette = {
                "A": QtGui.QColor(255, 214, 225),
                "C": QtGui.QColor(205, 232, 255),
                "G": QtGui.QColor(255, 235, 184),
                "U": QtGui.QColor(196, 243, 229),
                "T": QtGui.QColor(225, 211, 255),
                "I": QtGui.QColor(235, 211, 255),
            }
            return palette.get(base, QtGui.QColor(229, 235, 242))

        def paint(self, painter, option, widget=None):
            radius = float(self.RADIUS)
            gel = self.viewer.gel_style_enabled()
            saved = False
            try:
                painter.save()
                saved = True
                painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
                selected = bool(self.isSelected())
                hovered = bool(getattr(self, "_v4_hover", False))
                pressed = bool(getattr(self, "_v4_pressed", False))
                base = QtGui.QColor(self._base_fill())

                shadow = QtGui.QColor(2, 8, 23, 105 if gel else 38)
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QBrush(shadow))
                painter.drawEllipse(
                    QtCore.QRectF(
                        -radius + 2.2, -radius + 3.4,
                        2.0 * radius, 2.0 * radius,
                    )
                )

                if selected or hovered:
                    aura = QtGui.QColor(
                        55, 220, 255, 120 if selected else 55
                    )
                    aura_radius = radius + (5.2 if selected else 3.4)
                    painter.setPen(QtCore.Qt.NoPen)
                    painter.setBrush(QtGui.QBrush(aura))
                    painter.drawEllipse(
                        QtCore.QRectF(
                            -aura_radius, -aura_radius,
                            2.0 * aura_radius, 2.0 * aura_radius,
                        )
                    )

                if gel:
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
                    fill = QtGui.QBrush(base)
                    border = QtGui.QColor(55, 72, 96, 210)

                if selected:
                    border = QtGui.QColor(44, 215, 255, 255)
                elif hovered:
                    border = QtGui.QColor(105, 225, 255, 245)
                if pressed:
                    border = QtGui.QColor(255, 255, 255, 250)
                pen = QtGui.QPen(border)
                pen.setWidthF(2.35 if selected else (1.75 if hovered else 1.2))
                painter.setPen(pen)
                painter.setBrush(fill)
                painter.drawEllipse(
                    QtCore.QRectF(-radius, -radius, 2.0 * radius, 2.0 * radius)
                )

                if gel:
                    painter.setPen(QtCore.Qt.NoPen)
                    painter.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 145)))
                    painter.drawEllipse(
                        QtCore.QRectF(
                            -radius * 0.58, -radius * 0.68,
                            radius * 0.92, radius * 0.45,
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

        # Restore V4's spring machinery; it already checks jelly_cb at runtime.
        def _v4_set_target_scale(self, target, kick=0.0):
            return _DSSR_PY2D_NODE_BEFORE_V5._v4_set_target_scale(
                self, target, kick
            )

        def _v4_advance_visual(self):
            return _DSSR_PY2D_NODE_BEFORE_V5._v4_advance_visual(self)

        def hoverEnterEvent(self, event):
            return _DSSR_PY2D_NODE_BEFORE_V5.hoverEnterEvent(self, event)

        def hoverLeaveEvent(self, event):
            return _DSSR_PY2D_NODE_BEFORE_V5.hoverLeaveEvent(self, event)

        def itemChange(self, change, value):
            return _DSSR_PY2D_NODE_BEFORE_V5.itemChange(self, change, value)


    Dssr2DGraphicsView = Dssr2DGraphicsViewV6
    Dssr2DNodeItem = Dssr2DNodeItemV6


    class DssrPython2DDialogV6(_DSSR_PY2D_DIALOG_BEFORE_V6):
        """Premium RNA editor with switchable gel visuals and 2D/3D brushing."""

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
            self._v6_sync_timer = None
            self._v6_reverse_timer = None
            self._v6_sync_pending = False
            self._v6_sync_from_pymol = False
            self._v6_last_pymol_signature = None
            self._v6_last_highlight_signature = None
            self._v6_closed = False
            _DSSR_PY2D_DIALOG_BEFORE_V6.__init__(
                self,
                model=model,
                pymol_selection=pymol_selection,
                algorithm=algorithm,
                number_every=number_every,
                show_tertiary=show_tertiary,
                parent=parent,
            )
            self._v6_sync_timer = QtCore.QTimer(self)
            self._v6_sync_timer.setSingleShot(True)
            self._v6_sync_timer.setInterval(55)
            self._v6_sync_timer.timeout.connect(self._v6_flush_live_sync)
            self._v6_reverse_timer = QtCore.QTimer(self)
            self._v6_reverse_timer.setInterval(450)
            self._v6_reverse_timer.timeout.connect(self._v6_pull_pymol_selection)
            self._v6_reverse_timer.start()
            try:
                self.base_colors_cb.blockSignals(True)
                self.base_colors_cb.setChecked(False)
                self.base_colors_cb.blockSignals(False)
                self.base_colors = False
            except Exception:
                pass
            self.setWindowTitle("RNA 2D studio — %s" % self.model.title)
            self._v6_apply_appearance()
            self._v6_update_view_cursor()
            self._v4_ensure_animation()
            self._v3_update_editor_status("premium gel + 3D brush ready")

        def _v4_install_gel_controls(self):
            group = QtWidgets.QGroupBox("2D studio · appearance and interaction")
            outer = QtWidgets.QVBoxLayout(group)
            outer.setContentsMargins(10, 11, 10, 9)
            outer.setSpacing(7)

            top = QtWidgets.QHBoxLayout()
            top.addWidget(QtWidgets.QLabel("tool"))
            self.interaction_combo = QtWidgets.QComboBox()
            self.interaction_combo.addItem("edit / pan  [P]", "edit")
            self.interaction_combo.addItem("brush select  [B]", "brush")
            self.interaction_combo.setCurrentIndex(0)
            self.interaction_combo.setToolTip(
                "Brush: sweep across nucleotides. Edit/pan: drag bases or the canvas."
            )
            self.interaction_combo.currentIndexChanged.connect(
                lambda *_args: self._v6_interaction_changed()
            )
            top.addWidget(self.interaction_combo)

            top.addWidget(QtWidgets.QLabel("brush"))
            self.brush_radius_spin = QtWidgets.QSpinBox()
            self.brush_radius_spin.setRange(8, 100)
            self.brush_radius_spin.setValue(32)
            self.brush_radius_spin.setSuffix(" px")
            self.brush_radius_spin.setToolTip("Freehand brush radius on screen")
            top.addWidget(self.brush_radius_spin)

            self.gel_style_cb = QtWidgets.QCheckBox("gel mode")
            self.gel_style_cb.setChecked(True)
            self.gel_style_cb.setToolTip(
                "Toggle both the glass-like rendering and elastic motion"
            )
            self.gel_style_cb.toggled.connect(
                lambda checked: self._v6_gel_mode_toggled(checked)
            )
            top.addWidget(self.gel_style_cb)

            # Compatibility control consumed by V4's spring engine.  V6 uses
            # the single visible gel-mode switch above as the master control.
            self.jelly_cb = QtWidgets.QCheckBox()
            self.jelly_cb.setChecked(True)
            self.jelly_cb.setVisible(False)

            self.live_3d_cb = QtWidgets.QCheckBox("live 3D glow")
            self.live_3d_cb.setChecked(True)
            self.live_3d_cb.setToolTip(
                "Show selected 2D bases as a separate cyan object in PyMOL"
            )
            self.live_3d_cb.toggled.connect(
                lambda *_args: self._v3_sync_pymol_selection()
            )
            top.addWidget(self.live_3d_cb)

            self.reverse_3d_cb = QtWidgets.QCheckBox("3D → 2D sync")
            self.reverse_3d_cb.setChecked(True)
            self.reverse_3d_cb.setToolTip(
                "Mirror PyMOL's current 'sele' residues back into this 2D view"
            )
            self.reverse_3d_cb.toggled.connect(
                lambda checked: self._v6_reverse_sync_toggled(checked)
            )
            top.addWidget(self.reverse_3d_cb)

            self.zoom_3d_cb = QtWidgets.QCheckBox("zoom after brush")
            self.zoom_3d_cb.setChecked(False)
            top.addWidget(self.zoom_3d_cb)

            self.reload_2d_btn = QtWidgets.QPushButton("reload DSSR / state")
            self.reload_2d_btn.setToolTip(
                "Re-run DSSR for the object/state selected in the main GUI and reopen 2D"
            )
            self.reload_2d_btn.clicked.connect(self._v6_reload_from_pymol)
            top.addWidget(self.reload_2d_btn)
            top.addStretch(1)
            outer.addLayout(top)

            bottom = QtWidgets.QHBoxLayout()
            bottom.addWidget(QtWidgets.QLabel("drag mode"))
            self.drag_mode_combo = QtWidgets.QComboBox()
            for text, value in (
                ("soft neighborhood", "soft"),
                ("single base", "base"),
                ("current selection", "selection"),
                ("base pair", "pair"),
                ("loop / unpaired", "loop"),
                ("stem", "stem"),
                ("whole branch", "branch"),
            ):
                self.drag_mode_combo.addItem(text, value)
            self.drag_mode_combo.setCurrentIndex(0)
            self.drag_mode_combo.currentIndexChanged.connect(
                lambda *_args: self._v3_update_editor_status("drag mode changed")
            )
            bottom.addWidget(self.drag_mode_combo)

            self.follow_spin = QtWidgets.QDoubleSpinBox()
            self.follow_spin.setRange(0.10, 0.90)
            self.follow_spin.setSingleStep(0.05)
            self.follow_spin.setValue(0.62)
            self.follow_spin.setToolTip(
                "How strongly neighboring bases follow a dragged base in gel mode"
            )
            bottom.addWidget(QtWidgets.QLabel("elasticity"))
            bottom.addWidget(self.follow_spin)

            self.fit_selected_btn = QtWidgets.QPushButton("fit selected")
            self.fit_selected_btn.clicked.connect(self.fit_selected)
            bottom.addWidget(self.fit_selected_btn)
            bottom.addWidget(QtWidgets.QLabel("go to nt"))
            self.goto_spin = QtWidgets.QSpinBox()
            self.goto_spin.setRange(1, max(1, len(self.model.nts)))
            self.goto_spin.setValue(1)
            bottom.addWidget(self.goto_spin)
            self.goto_btn = QtWidgets.QPushButton("go")
            self.goto_btn.clicked.connect(self.goto_nucleotide)
            bottom.addWidget(self.goto_btn)

            hint = QtWidgets.QLabel(
                "Edit: drag bases normally · B=brush · Shift+brush=add · "
                "Alt+brush=erase · P=return to edit · right/middle drag=pan"
            )
            hint.setWordWrap(True)
            hint.setObjectName("studioHint")
            bottom.addWidget(hint, 1)
            outer.addLayout(bottom)

            root = self.layout()
            if root is not None:
                try:
                    root.insertWidget(2, group)
                except Exception:
                    root.addWidget(group)
            self.gel_group = group

        def interaction_tool(self):
            try:
                return str(self.interaction_combo.currentData() or "edit")
            except Exception:
                return "edit"

        def set_interaction_tool(self, tool):
            wanted = str(tool or "edit").lower()
            try:
                for index in range(self.interaction_combo.count()):
                    if str(self.interaction_combo.itemData(index)) == wanted:
                        self.interaction_combo.setCurrentIndex(index)
                        return
            except Exception:
                pass

        def _v6_interaction_changed(self):
            self._v6_update_view_cursor()
            self._v3_update_editor_status("tool=%s" % self.interaction_tool())

        def _v6_update_view_cursor(self):
            try:
                if self.interaction_tool() == "brush":
                    self.view.setCursor(QtCore.Qt.CrossCursor)
                else:
                    self.view.setCursor(QtCore.Qt.OpenHandCursor)
            except Exception:
                pass

        def gel_style_enabled(self):
            try:
                return bool(self.gel_style_cb.isChecked())
            except Exception:
                return True

        def _v6_gel_mode_toggled(self, checked):
            try:
                self.jelly_cb.setChecked(bool(checked))
            except Exception:
                pass
            self._v6_apply_appearance()
            self._v4_ensure_animation()
            self._v3_update_editor_status(
                "gel mode on" if checked else "gel mode off"
            )

        def _v4_apply_visual_theme(self):
            self._v6_apply_appearance()

        def _v6_apply_appearance(self):
            style = """
                QDialog { background: #ffffff; color: #263746; }
                QLabel, QCheckBox { color: #263746; }
                QGroupBox {
                    color: #263746; border: 1px solid #d4dde5;
                    border-radius: 10px; margin-top: 8px; padding-top: 9px;
                    background: #ffffff; font-weight: 600;
                }
                QGroupBox::title { subcontrol-origin: margin; left: 10px;
                    padding: 0 6px; color: #24657d; }
                QPushButton, QComboBox, QSpinBox, QDoubleSpinBox {
                    color: #263746; background: #ffffff;
                    border: 1px solid #bdcbd6; border-radius: 7px;
                    padding: 4px 8px; min-height: 22px;
                }
                QPushButton:hover, QComboBox:hover {
                    background: #eef9fd; border-color: #52b9da;
                }
                QPushButton:pressed { background: #dceff6; }
                QComboBox QAbstractItemView {
                    color: #263746; background: #ffffff;
                    selection-color: #123d50;
                    selection-background-color: #d6f0fa;
                }
                QCheckBox::indicator { width: 15px; height: 15px; }
                QLabel#studioHint { color: #607682; font-weight: 400; }
                QToolTip { color: #263746; background: #ffffff;
                    border: 1px solid #52b9da; }
            """
            canvas_color = QtGui.QColor(255, 255, 255)
            try:
                self.setStyleSheet(style)
                self.view.setBackgroundBrush(QtGui.QBrush(canvas_color))
            except Exception:
                pass
            self._v6_refresh_scene_style()

        def _v6_refresh_scene_style(self):
            gel = self.gel_style_enabled()
            for node in list(getattr(self, "nodes", [])):
                try:
                    node.base_text_item.setBrush(
                        QtGui.QBrush(
                            QtGui.QColor(16, 42, 59)
                            if gel else QtGui.QColor(31, 50, 68)
                        )
                    )
                    for child in node.childItems():
                        if child is node.base_text_item:
                            continue
                        if isinstance(child, QtWidgets.QGraphicsSimpleTextItem):
                            child.setBrush(
                                QtGui.QBrush(
                                    QtGui.QColor(72, 94, 112)
                                )
                            )
                    node.update()
                except Exception:
                    pass
            try:
                self.scene.update()
            except Exception:
                pass

        def _v4_ensure_animation(self):
            return _DSSR_PY2D_DIALOG_BEFORE_V5._v4_ensure_animation(self)

        def _v4_animation_tick(self):
            return _DSSR_PY2D_DIALOG_BEFORE_V5._v4_animation_tick(self)

        def redraw(self, *_args):
            _DSSR_PY2D_DIALOG_BEFORE_V6.redraw(self, *_args)
            self._v6_apply_appearance()
            self._v4_ensure_animation()

        def _v3_selection_changed(self):
            if getattr(self, "_v3_rebuilding", False):
                return
            self._v6_schedule_live_sync()
            self._v3_update_editor_status("selection changed")

        def _v6_after_brush_selection(self, final=False):
            selected = sorted(
                node.nt_index for node in self.nodes if node.isSelected()
            )
            self._v3_active_index = selected[0] if selected else None
            if final:
                self._v6_flush_live_sync(final=True)
                self._v3_update_editor_status("brush selection mapped to 3D")
            else:
                self._v6_schedule_live_sync()
                self._v3_update_editor_status("brushing")

        def _v6_schedule_live_sync(self):
            if getattr(self, "_v6_closed", False):
                return
            self._v6_sync_pending = True
            timer = getattr(self, "_v6_sync_timer", None)
            if timer is None:
                self._v6_flush_live_sync()
                return
            try:
                timer.start()
            except Exception:
                self._v6_flush_live_sync()

        def _v6_flush_live_sync(self, final=False):
            if getattr(self, "_v6_closed", False):
                return
            self._v6_sync_pending = False
            self._v3_sync_pymol_selection()
            if final:
                try:
                    if self.zoom_3d_cb.isChecked() and cmd.count_atoms(
                        self.HIGHLIGHT_OBJECT
                    ):
                        cmd.zoom(self.HIGHLIGHT_OBJECT, buffer=4.0)
                except Exception:
                    pass

        def _v6_node_residue_signature(self):
            residues = set()
            for node in getattr(self, "nodes", []):
                if not node.isSelected():
                    continue
                nt = self.model.nts[node.nt_index]
                chain = str(nt.get("chain", "")).strip()
                resi = str(nt.get("resi", "")).strip()
                if resi:
                    residues.add((chain, resi))
            return tuple(sorted(residues))

        def _v6_reverse_sync_toggled(self, checked):
            self._v6_last_pymol_signature = None
            if checked:
                self._v6_pull_pymol_selection()
                self._v3_update_editor_status("bidirectional sync on")
            else:
                self._v3_update_editor_status("3D-to-2D sync off")

        def _v6_pull_pymol_selection(self):
            if (
                getattr(self, "_v6_closed", False)
                or getattr(self, "_v6_sync_pending", False)
                or getattr(self, "_v6_sync_from_pymol", False)
            ):
                return
            try:
                if not self.isVisible():
                    return
            except Exception:
                pass
            try:
                if not self.reverse_3d_cb.isChecked():
                    return
            except Exception:
                return

            residues = set()
            try:
                scoped = "((%s) and sele)" % self.pymol_selection
                if int(cmd.count_atoms(scoped)) > 0:
                    cmd.iterate(
                        scoped,
                        "_dssr_residues.add((chain, resi))",
                        space={"_dssr_residues": residues},
                    )
            except Exception:
                residues = set()

            signature = tuple(sorted(residues))
            if signature == self._v6_last_pymol_signature:
                return
            self._v6_last_pymol_signature = signature

            wanted = set()
            for index, nt in enumerate(getattr(self.model, "nts", [])):
                key = (
                    str(nt.get("chain", "")).strip(),
                    str(nt.get("resi", "")).strip(),
                )
                if key[1] and key in residues:
                    wanted.add(index)

            current = {
                node.nt_index
                for node in getattr(self, "nodes", [])
                if node.isSelected()
            }
            if wanted == current:
                return

            self._v6_sync_from_pymol = True
            self._v3_rebuilding = True
            try:
                for node in getattr(self, "nodes", []):
                    node.setSelected(node.nt_index in wanted)
                self._v3_active_index = min(wanted) if wanted else None
            finally:
                self._v3_rebuilding = False
                self._v6_sync_from_pymol = False
            self._v6_update_pymol_highlight()
            self._v3_update_editor_status("3D selection mirrored to 2D")

        def _v6_reload_from_pymol(self):
            parent = self.parent()
            if parent is None or not hasattr(parent, "_get_dssr_context"):
                self._v3_update_editor_status(
                    "reload is available when opened from dssr_gui"
                )
                return
            try:
                parent._invalidate_dssr_cache()
                replacement = DssrPython2DIntegration.gui_open(parent)
                if replacement is not None:
                    self._v3_update_editor_status("reloaded DSSR context")
                    QtCore.QTimer.singleShot(0, self.close)
            except Exception as error:
                self._v3_update_editor_status("reload error: %s" % str(error))

        def _v3_sync_pymol_selection(self):
            _DSSR_PY2D_DIALOG_BEFORE_V5._v3_sync_pymol_selection(self)
            self._v6_last_pymol_signature = self._v6_node_residue_signature()
            self._v6_update_pymol_highlight()

        def _v6_update_pymol_highlight(self):
            name = self.HIGHLIGHT_OBJECT
            try:
                enabled = bool(self.live_3d_cb.isChecked())
            except Exception:
                enabled = False
            signature = self._v6_node_residue_signature()
            should_show = (
                not getattr(self, "_v6_closed", False)
                and enabled
                and bool(signature)
            )
            if should_show:
                try:
                    object_exists = name in cmd.get_names("all")
                except Exception:
                    object_exists = False
                if (
                    signature == self._v6_last_highlight_signature
                    and object_exists
                ):
                    return
            else:
                object_exists = False

            try:
                cmd.delete(name)
                _DSSR_BLOCK_OBJECTS.discard(name)
            except Exception:
                pass
            self._v6_last_highlight_signature = None
            if not should_show:
                return
            try:
                if int(cmd.count_atoms("sele")) <= 0:
                    return
                try:
                    source_state = max(1, int(cmd.get_state()))
                except Exception:
                    source_state = 1
                cmd.create(
                    name,
                    "byres (sele)",
                    source_state=source_state,
                    target_state=1,
                    zoom=0,
                    quiet=1,
                )
                _DSSR_BLOCK_OBJECTS.add(name)
                self._v6_last_highlight_signature = signature
                cmd.hide("everything", name)
                cmd.show("sticks", name)
                cmd.show("spheres", "(%s) and name P" % name)
                cmd.set_color("dssr_2d_glow", [0.08, 0.86, 1.00])
                cmd.color("dssr_2d_glow", name)
                cmd.set("stick_radius", 0.23, name)
                cmd.set("stick_transparency", 0.08, name)
                cmd.set("sphere_scale", 0.34, name)
                cmd.enable(name)
                cmd.refresh()
            except Exception as error:
                try:
                    self.status_label.setText(
                        "3D highlight error: %s" % str(error)
                    )
                except Exception:
                    pass

        def _v3_update_editor_status(self, action=""):
            selected = sum(
                1 for node in getattr(self, "nodes", []) if node.isSelected()
            )
            variant = str(
                getattr(self.model, "_dssr2d_layout_variant", self.algorithm)
            )
            text = (
                self.model.summary()
                + " | %s" % variant
                + " | selected=%d" % selected
                + " | tool=%s" % self.interaction_tool()
                + " | B=brush, P=edit/pan, Shift=add, Alt=erase"
            )
            if action:
                text += " | %s" % action
            try:
                self.status_label.setText(text)
            except Exception:
                pass

        def closeEvent(self, event):
            self._v6_closed = True
            self._v6_sync_pending = False
            try:
                if self._v6_sync_timer is not None:
                    self._v6_sync_timer.stop()
                if self._v6_reverse_timer is not None:
                    self._v6_reverse_timer.stop()
                if hasattr(self, "_v4_timer"):
                    self._v4_timer.stop()
            except Exception:
                pass
            try:
                cmd.delete(self.HIGHLIGHT_OBJECT)
                _DSSR_BLOCK_OBJECTS.discard(self.HIGHLIGHT_OBJECT)
            except Exception:
                pass
            try:
                QtWidgets.QDialog.closeEvent(self, event)
            except Exception:
                event.accept()


    DssrPython2DDialog = DssrPython2DDialogV6
    DssrPython2DDialog._premium_editor_v6_installed = True


# Refresh the labels in the main dssr_gui window without changing its flow.
_DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V6 = (
    DssrPython2DIntegration.install_gui_controls
)


def _dssr2d_v6_install_gui_controls(dialog):
    _DSSR_PY2D_INSTALL_CONTROLS_BEFORE_V6(dialog)
    try:
        dialog.py2d_group.setTitle(
            "RNA 2D studio (gel + free brush + bidirectional 3D sync)"
        )
        dialog.py2d_open_btn.setToolTip(
            "Open the pure-Python RNA 2D studio. Toggle gel styling, sweep "
            "across bases with the brush, and synchronize selections both ways."
        )
    except Exception:
        pass


DssrPython2DIntegration.install_gui_controls = staticmethod(
    _dssr2d_v6_install_gui_controls
)
DssrPython2DIntegration._premium_editor_v6_installed = True

try:
    print(
        "Loaded DSSR RNA 2D premium editor %s: switchable gel, free brush, "
        "bidirectional 2D/3D selection"
        % __DSSR_PY2D_PREMIUM_EDITOR_VERSION__
    )
except Exception:
    pass
