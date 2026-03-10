# dssr_select.py
# DSSR-based selection of RNA structural features in PyMOL
#
# (c) 2026 Eric Chen, Bener Dulger, and Xiang-Jun Lu
#
# This project was initiated and coordinated by Xiang-Jun Lu.
#
# CONTRIBUTIONS:
# - Eric Chen: Developed the Qt-based GUI, incorporated dssr_block functionality,
#   enhanced feature parsing, and performed final code consolidation.
# - Bener Dulger: Initial implementation of structural feature selection and JSON parsing.
# - Thomas Holder: Original 'dssr_block' logic (c) Schrodinger LLC.
#
# LICENSE: BSD 2-Clause
#
# This plugin incorporates code from 'dssr_block'.
# Redistributions must retain the original copyright notice and this license.

from pymol import cmd, CmdException

__DSSR_PLUGIN_VERSION__ = 'v1.0.0'

_hex_color_cache = {}
selected_features = []
_DSSR_BLOCK_OBJECTS = set()

FEATURE_MAP = {
    'pairs': 'pairs',
    'hbonds': 'hbonds',
    'stems': 'stems',
    'helices': 'helices',
    'stacks': 'stacks',
    'nonstack': 'nonStack',
    'coaxstacks': 'coaxStacks',
    'atom2bases': 'atom2bases',
    'aminors': 'Aminors',
    'splayunits': 'splayUnits',
    'hairpins': 'hairpins',
    'bulges': 'bulges',
    'iloops': 'iloops',
    'internal': 'iloops',
    'junctions': 'junctions',
    'sssegments': 'ssSegments',
    'ssSegments': 'ssSegments',
    'multiplets': 'multiplets',
    'nts': 'nts',
    'pseudoknot': 'pseudoknot',
}

FEATURE_ORDER = [
    'pairs', 'hbonds', 'stems', 'helices',
    'stacks', 'nonstack', 'hairpins', 'bulges',
    'iloops', 'junctions', 'sssegments', 'multiplets',
    'coaxstacks', 'atom2bases', 'aminors', 'splayunits',
    'nts', 'pseudoknot'
]


def unquote(s):
    s = str(s)
    if not s:
        return s
    if s.rstrip()[-1:] not in ('"', "'"):
        return s
    return cmd.safe_eval(s)


def _safe_tail(s, n=500):
    try:
        s = str(s)
    except Exception:
        return ''
    if len(s) <= n:
        return s
    return s[-n:]


def run_dssr_json(pdb_path, exe):
    import subprocess, json
    args = [exe, '--json', '-i=' + pdb_path]

    try:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        rc = p.returncode
    except OSError:
        raise CmdException('Cannot execute exe="%s"' % exe)

    try:
        out_txt = out.decode('utf-8', errors='replace') if out else ''
    except Exception:
        out_txt = str(out)
    try:
        err_txt = err.decode('utf-8', errors='replace') if err else ''
    except Exception:
        err_txt = str(err)

    if rc != 0:
        raise CmdException('DSSR failed (rc=%s). stderr tail: %s' % (str(rc), _safe_tail(err_txt)))

    if not out_txt.strip():
        raise CmdException('DSSR returned empty stdout (expected JSON). stderr tail: %s' % _safe_tail(err_txt))

    try:
        return json.loads(out_txt)
    except Exception:
        s = out_txt
        i = s.find('{')
        j = s.rfind('}')
        if i >= 0 and j > i:
            try:
                return json.loads(s[i:j + 1])
            except Exception as e2:
                raise CmdException('Failed to parse DSSR JSON. stdout head: %s | stderr tail: %s | err: %s' %
                                   (s[:120].replace('\n', ' '), _safe_tail(err_txt), str(e2)))
        raise CmdException('Failed to parse DSSR JSON (no JSON object found). stdout head: %s | stderr tail: %s' %
                           (s[:120].replace('\n', ' '), _safe_tail(err_txt)))


def _hex_to_rgb01(h):
    import re
    s = str(h).strip()

    if s.startswith('"') and s.endswith('"'):
        s = s[1:-1].strip()
    if s.startswith("'") and s.endswith("'"):
        s = s[1:-1].strip()

    if s.lower().startswith('0x'):
        s = s[2:]
    if s.startswith('#'):
        s = s[1:]

    s = s.strip()
    if re.fullmatch(r'[0-9a-fA-F]{3}', s):
        s = ''.join([c * 2 for c in s])

    if not re.fullmatch(r'[0-9a-fA-F]{6}', s):
        raise CmdException('Invalid hex color "%s". Use FF00AA or 0xFF00AA (or "#FF00AA" quoted).' % h)

    r = int(s[0:2], 16) / 255.0
    g = int(s[2:4], 16) / 255.0
    b = int(s[4:6], 16) / 255.0
    return s.lower(), [r, g, b]


def _resolve_color_spec(color_spec):
    if color_spec is None:
        return None
    s = str(color_spec).strip()
    if not s:
        return None
    if s.lower() in ('auto', 'default'):
        return None

    is_hexish = (
        s.startswith('#') or
        s.lower().startswith('0x') or
        (len(s) in (3, 6) and all(c in '0123456789abcdefABCDEF' for c in s))
    )
    if is_hexish:
        hex6, rgb = _hex_to_rgb01(s)
        if hex6 in _hex_color_cache:
            return _hex_color_cache[hex6]
        cname = 'dssr_hex_%s' % hex6
        cmd.set_color(cname, rgb)
        _hex_color_cache[hex6] = cname
        return cname

    name = s.lower()
    try:
        idx = cmd.get_color_index(name)
        if idx < 0:
            raise Exception('unknown')
    except Exception:
        raise CmdException('Unknown color "%s". Use a PyMOL color name (e.g., blue) or hex like FF00AA / 0xFF00AA.' % s)
    return name


def parse_nt_id(nt_id):
    import re
    s = str(nt_id).strip()

    if '.' in s and not s.startswith('///'):
        try:
            chain, res = s.split('.', 1)
        except ValueError:
            raise CmdException('Unexpected nt_id format: "%s"' % s)
        m = re.search(r'(-?\d+)$', res)
        if not m:
            raise CmdException('Could not extract residue number from "%s"' % res)
        return chain, m.group(1)

    if '/' in s:
        parts = s.strip('/').split('/')
        if len(parts) < 2:
            raise CmdException('Unexpected nt_id format: "%s"' % s)
        chain = parts[-2]
        res = parts[-1]
        m = re.search(r'(-?\d+)$', res)
        if not m:
            raise CmdException('Could not extract residue number from "%s"' % res)
        return chain, m.group(1)

    raise CmdException('Unexpected nt_id format: "%s"' % s)


def parse_atom_id(atom_id):
    try:
        _, rest = atom_id.split('@', 1)
    except ValueError:
        raise CmdException('Unexpected atom_id format: "%s"' % atom_id)
    return parse_nt_id(rest)


def parse_hbond_atom(atom_id):
    import re
    s = str(atom_id).strip()

    if '@' in s:
        atom_name, nt = s.split('@', 1)
        chain, resi = parse_nt_id(nt)
        return chain, resi, atom_name

    if '/' in s:
        parts = s.strip('/').split('/')
        if len(parts) < 3:
            raise CmdException('Unexpected hbond atom format: "%s"' % s)
        chain = parts[-3]
        res = parts[-2]
        atom_name = parts[-1]
        m = re.search(r'(-?\d+)$', res)
        if not m:
            raise CmdException('Could not extract residue number from "%s"' % res)
        return chain, m.group(1), atom_name

    raise CmdException('Unexpected hbond atom format: "%s"' % s)


def parse_a2b_atom(atom_id):
    import re
    s = str(atom_id).strip()

    if '@' in s:
        atom_name, nt = s.split('@', 1)
        chain, resi = parse_nt_id(nt)
        return chain, resi, atom_name

    if '/' in s:
        parts = s.strip('/').split('/')
        if len(parts) < 3:
            raise CmdException('Unexpected atom format: "%s"' % s)
        chain = parts[-3]
        res = parts[-2]
        atom_name = parts[-1]
        m = re.search(r'(-?\d+)$', res)
        if not m:
            raise CmdException('Could not extract residue number from "%s"' % res)
        return chain, m.group(1), atom_name

    raise CmdException('Unexpected atom format: "%s"' % s)


def parse_dotbracket_pseudoknots(dotbracket):
    stack_map = {'(': ')', '[': ']', '{': '}', '<': '>'}
    close_to_open = {v: k for k, v in stack_map.items()}
    stacks = {}
    layers = {}
    layer_assignment = {'()': 0, '[]': 1, '{}': 2, '<>': 3}
    next_layer = 4

    for idx, char in enumerate(dotbracket):
        if char == '.':
            continue

        if char in stack_map:
            bracket_type = char + stack_map[char]
            stacks.setdefault(bracket_type, []).append(idx)

        elif char in close_to_open:
            open_char = close_to_open[char]
            bracket_type = open_char + char
            if bracket_type not in stacks or not stacks[bracket_type]:
                continue
            open_idx = stacks[bracket_type].pop()
            if bracket_type not in layer_assignment:
                layer_assignment[bracket_type] = next_layer
                next_layer += 1
            layer = layer_assignment[bracket_type]
            layers.setdefault(layer, []).append((open_idx, idx))

        elif char.isalpha():
            stacks.setdefault(char, [])
            if not stacks[char]:
                stacks[char].append(idx)
            else:
                open_idx = stacks[char].pop()
                if char not in layer_assignment:
                    layer_assignment[char] = next_layer
                    next_layer += 1
                layer = layer_assignment[char]
                layers.setdefault(layer, []).append((open_idx, idx))

    return layers


def _atom_sel(chain, resi, atom_name):
    atom_name = str(atom_name).replace('"', '\\"')
    return '(chain %s and resi %s and name "%s")' % (chain, resi, atom_name)


def build_selection_from_layer(layer_pairs, nts_list):
    residues = set()
    for open_idx, close_idx in layer_pairs:
        if open_idx < len(nts_list):
            nt1_id = nts_list[open_idx].get('nt_id')
            if nt1_id:
                residues.add(parse_nt_id(nt1_id))
        if close_idx < len(nts_list):
            nt2_id = nts_list[close_idx].get('nt_id')
            if nt2_id:
                residues.add(parse_nt_id(nt2_id))

    if not residues:
        return None

    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


def build_selection_from_pair(pair_entry):
    nt1 = pair_entry.get('nt1')
    nt2 = pair_entry.get('nt2')
    if not nt1 or not nt2:
        raise CmdException('Pair entry missing nt1 or nt2')
    residues = {parse_nt_id(nt1), parse_nt_id(nt2)}
    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


def build_selection_from_hbond(hb_entry, hbonds_mode='residue'):
    atom1_id = hb_entry.get('atom1_id')
    atom2_id = hb_entry.get('atom2_id')
    if not atom1_id or not atom2_id:
        raise CmdException('H-bond entry missing atom1_id or atom2_id')

    mode = str(hbonds_mode).strip().lower()
    if mode in ('atom', 'distance'):
        c1, r1, a1 = parse_hbond_atom(atom1_id)
        c2, r2, a2 = parse_hbond_atom(atom2_id)
        s1 = _atom_sel(c1, r1, a1)
        s2 = _atom_sel(c2, r2, a2)
        return ('%s or %s' % (s1, s2), s1, s2)

    residues = {parse_atom_id(atom1_id), parse_atom_id(atom2_id)}
    sel = ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)
    return (sel, None, None)


def build_selection_from_nts_list(nts_list):
    if not nts_list:
        raise CmdException('Empty nucleotide list')
    residues = {parse_nt_id(nt) for nt in nts_list}
    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


def build_selection_from_stem(stem_entry):
    pairs = stem_entry.get('pairs', [])
    if not pairs:
        raise CmdException('Stem has no pairs')

    residues = set()
    for p in pairs:
        if p.get('nt1'):
            residues.add(parse_nt_id(p['nt1']))
        if p.get('nt2'):
            residues.add(parse_nt_id(p['nt2']))

    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


def build_selection_from_hairpin(hairpin_entry):
    nts_long = hairpin_entry.get('nts_long')
    if not nts_long:
        raise CmdException('Hairpin missing nts_long field')
    nts_list = [nt.strip() for nt in nts_long.split(',') if nt.strip()]
    return build_selection_from_nts_list(nts_list)


def build_selection_from_coaxstack(coax_entry, stems_list):
    stem_indices = coax_entry.get('stem_indices', [])
    if not stem_indices:
        raise CmdException('coaxStacks entry missing stem_indices')

    residues = set()
    for si in stem_indices:
        try:
            idx = int(si)
        except Exception:
            continue
        if idx < 1 or idx > len(stems_list):
            continue
        stem_entry = stems_list[idx - 1]
        pairs = stem_entry.get('pairs', [])
        for p in pairs:
            if p.get('nt1'):
                residues.add(parse_nt_id(p['nt1']))
            if p.get('nt2'):
                residues.add(parse_nt_id(p['nt2']))

    if not residues:
        raise CmdException('Could not build selection for coaxStacks entry')

    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


def build_selection_from_atom2base(a2b_entry):
    atom = a2b_entry.get('atom')
    nt = a2b_entry.get('nt')

    clauses = []

    if nt:
        c_nt, r_nt = parse_nt_id(nt)
        clauses.append('(chain %s and resi %s)' % (c_nt, r_nt))

    if atom:
        c_a, r_a, atom_name = parse_a2b_atom(atom)
        atom_name = str(atom_name).replace('"', '\\"')
        clauses.append('(chain %s and resi %s and name "%s")' % (c_a, r_a, atom_name))

    if not clauses:
        raise CmdException('atom2bases entry missing atom and nt')

    return ' or '.join(clauses)


def build_selection_from_aminor(aminor_entry):
    desc_long = aminor_entry.get('desc_long', '')
    if not desc_long or 'vs' not in desc_long:
        raise CmdException('Aminors entry missing desc_long')

    left, right = desc_long.split('vs', 1)
    nts = []
    left = left.strip()
    right = right.strip()

    if left:
        nts.append(left)
    if right:
        for item in right.split(','):
            item = item.strip()
            if item:
                nts.append(item)

    if not nts:
        raise CmdException('Aminors entry has empty residues')

    return build_selection_from_nts_list(nts)


def _shorten_nts_long(nts_long, max_items=6):
    if not nts_long:
        return ''
    nts = [x.strip() for x in nts_long.split(',') if x.strip()]
    if len(nts) <= max_items:
        return ', '.join(nts)
    half = max_items // 2
    return '%s, ..., %s' % (', '.join(nts[:half]), ', '.join(nts[-half:]))


def _preview_entry(feature, entry, i):
    if feature == 'pairs':
        nt1 = entry.get('nt1', '?')
        nt2 = entry.get('nt2', '?')
        lw = entry.get('LW', entry.get('bp', ''))
        return '%d: %s - %s%s' % (i, nt1, nt2, (' (%s)' % lw) if lw else '')

    if feature == 'hbonds':
        return '%d: %s -- %s' % (i, entry.get('atom1_id', '?'), entry.get('atom2_id', '?'))

    if feature in ('stems', 'helices'):
        n = len(entry.get('pairs', [])) if isinstance(entry.get('pairs', []), list) else 0
        nm = entry.get('name', entry.get('index', ''))
        return '%d: %s (pairs=%d)' % (i, str(nm), n)

    if feature in ('stacks', 'nonstack', 'hairpins', 'bulges', 'iloops', 'internal', 'junctions',
                   'sssegments', 'ssSegments', 'multiplets', 'splayunits'):
        s = _shorten_nts_long(entry.get('nts_long', ''))
        return '%d: %s' % (i, s if s else '(missing nts_long)')

    if feature == 'coaxstacks':
        return '%d: helix=%s stems=%s' % (i, str(entry.get('helix_index', '')), str(entry.get('stem_indices', [])))

    if feature == 'atom2bases':
        t = entry.get('type', '')
        atom = entry.get('atom', '?')
        nt = entry.get('nt', '?')
        return '%d: %s atom=%s nt=%s' % (i, t if t else 'entry', atom, nt)

    if feature == 'aminors':
        ds = entry.get('desc_short', '')
        dl = entry.get('desc_long', '')
        return '%d: %s' % (i, ds if ds else dl)

    if feature == 'nts':
        return '%d: %s' % (i, entry.get('nt_id', '?'))

    return '%d: (no preview)' % i


def _extract_dotbracket(dssr_data):
    if 'dbn' not in dssr_data:
        raise CmdException('No dot-bracket notation found in DSSR output')
    dbn_data = dssr_data['dbn']

    if isinstance(dbn_data, dict):
        if isinstance(dbn_data.get('all_chains'), dict) and 'sstr' in dbn_data['all_chains']:
            return dbn_data['all_chains']['sstr']
        if 'sstr' in dbn_data:
            return dbn_data['sstr']
        for v in dbn_data.values():
            if isinstance(v, dict) and 'sstr' in v:
                return v['sstr']
        raise CmdException('Could not find sstr field in dbn data')

    return dbn_data


def _extract_chain_names(dssr_data):
    chains = set()
    try:
        nts = dssr_data.get('nts', [])
    except Exception:
        nts = []
    if isinstance(nts, list):
        for nt in nts:
            try:
                nt_id = nt.get('nt_id', '')
            except Exception:
                nt_id = ''
            if nt_id:
                try:
                    c, _ = parse_nt_id(nt_id)
                except Exception:
                    c = None
                if c:
                    chains.add(str(c))
    return sorted(chains)


def _count_pseudoknot_layers(dssr_data):
    try:
        dotbracket = _extract_dotbracket(dssr_data)
        layers = parse_dotbracket_pseudoknots(dotbracket)
        return len(layers) if layers else 0
    except Exception:
        return 0


def _format_rna_summary_text(dssr_data):
    chains = _extract_chain_names(dssr_data)
    chains_txt = ' '.join(chains) if chains else '(unknown)'

    pairs_n = len(dssr_data.get('pairs', [])) if isinstance(dssr_data.get('pairs', None), list) else 0
    hairpins_n = len(dssr_data.get('hairpins', [])) if isinstance(dssr_data.get('hairpins', None), list) else 0
    stems_n = len(dssr_data.get('stems', [])) if isinstance(dssr_data.get('stems', None), list) else 0
    bulges_n = len(dssr_data.get('bulges', [])) if isinstance(dssr_data.get('bulges', None), list) else 0
    junctions_n = len(dssr_data.get('junctions', [])) if isinstance(dssr_data.get('junctions', None), list) else 0
    pk_n = _count_pseudoknot_layers(dssr_data)
    aminors_n = len(dssr_data.get('Aminors', [])) if isinstance(dssr_data.get('Aminors', None), list) else 0
    stacks_n = len(dssr_data.get('stacks', [])) if isinstance(dssr_data.get('stacks', None), list) else 0

    lines = []
    lines.append('RNA Structure Summary')
    lines.append('---------------------')
    lines.append('Chains: %s' % chains_txt)
    lines.append('Base pairs: %d' % pairs_n)
    lines.append('Hairpins: %d' % hairpins_n)
    lines.append('Stems: %d' % stems_n)
    lines.append('Bulges: %d' % bulges_n)
    lines.append('Junctions: %d' % junctions_n)
    lines.append('Pseudoknots: %d' % pk_n)
    lines.append('A-minor interactions: %d' % aminors_n)
    lines.append('Stacking interactions: %d' % stacks_n)
    return '\n'.join(lines)


def _collect_residues_from_nts_long(nts_long):
    residues = set()
    if not nts_long:
        return residues
    nts_list = [nt.strip() for nt in str(nts_long).split(',') if nt.strip()]
    for nt in nts_list:
        try:
            residues.add(parse_nt_id(nt))
        except Exception:
            pass
    return residues


def _collect_residues_all(dssr_data, feature):
    feature = str(feature).lower().strip()
    residues = set()

    if feature == 'pairs':
        pairs = dssr_data.get('pairs', [])
        if isinstance(pairs, list):
            for p in pairs:
                try:
                    nt1 = p.get('nt1')
                    nt2 = p.get('nt2')
                except Exception:
                    nt1, nt2 = None, None
                if nt1:
                    try:
                        residues.add(parse_nt_id(nt1))
                    except Exception:
                        pass
                if nt2:
                    try:
                        residues.add(parse_nt_id(nt2))
                    except Exception:
                        pass
        return residues

    if feature == 'stems':
        stems = dssr_data.get('stems', [])
        if isinstance(stems, list):
            for st in stems:
                try:
                    pairs = st.get('pairs', [])
                except Exception:
                    pairs = []
                if not isinstance(pairs, list):
                    continue
                for p in pairs:
                    try:
                        nt1 = p.get('nt1')
                        nt2 = p.get('nt2')
                    except Exception:
                        nt1, nt2 = None, None
                    if nt1:
                        try:
                            residues.add(parse_nt_id(nt1))
                        except Exception:
                            pass
                    if nt2:
                        try:
                            residues.add(parse_nt_id(nt2))
                        except Exception:
                            pass
        return residues

    if feature in ('hairpins', 'bulges', 'junctions', 'stacks'):
        json_key = FEATURE_MAP.get(feature, feature)
        entries = dssr_data.get(json_key, [])
        if isinstance(entries, list):
            for e in entries:
                try:
                    nts_long = e.get('nts_long', '')
                except Exception:
                    nts_long = ''
                residues |= _collect_residues_from_nts_long(nts_long)
        return residues

    if feature == 'aminors':
        entries = dssr_data.get('Aminors', [])
        if isinstance(entries, list):
            for a in entries:
                try:
                    desc_long = a.get('desc_long', '')
                except Exception:
                    desc_long = ''
                if not desc_long or 'vs' not in desc_long:
                    continue
                try:
                    left, right = desc_long.split('vs', 1)
                except Exception:
                    continue
                nts = []
                left = left.strip()
                right = right.strip()
                if left:
                    nts.append(left)
                if right:
                    for item in right.split(','):
                        item = item.strip()
                        if item:
                            nts.append(item)
                for nt in nts:
                    try:
                        residues.add(parse_nt_id(nt))
                    except Exception:
                        pass
        return residues

    if feature == 'pseudoknot':
        try:
            dotbracket = _extract_dotbracket(dssr_data)
            nts_list = dssr_data.get('nts', None)
            if nts_list is None or not isinstance(nts_list, list):
                return residues
            layers = parse_dotbracket_pseudoknots(dotbracket) or {}
            for pairs in layers.values():
                for open_idx, close_idx in pairs:
                    if open_idx < len(nts_list):
                        nt1_id = nts_list[open_idx].get('nt_id')
                        if nt1_id:
                            try:
                                residues.add(parse_nt_id(nt1_id))
                            except Exception:
                                pass
                    if close_idx < len(nts_list):
                        nt2_id = nts_list[close_idx].get('nt_id')
                        if nt2_id:
                            try:
                                residues.add(parse_nt_id(nt2_id))
                            except Exception:
                                pass
        except Exception:
            pass
        return residues

    return residues


def _sort_resi_key(resi_str):
    try:
        return (0, int(str(resi_str)))
    except Exception:
        return (1, str(resi_str))


def _compact_sel_from_residues(residues):
    if not residues:
        return ''
    by_chain = {}
    for c, r in residues:
        c = str(c)
        r = str(r)
        by_chain.setdefault(c, set()).add(r)

    parts = []
    for c in sorted(by_chain.keys()):
        resis = sorted(by_chain[c], key=_sort_resi_key)
        resi_expr = '+'.join(resis)
        parts.append('(chain %s and resi %s)' % (c, r))
    return ' or '.join(parts)

def _build_residue_sel_from_dssr(dssr_data, feature, index):
    feature = str(feature).lower().strip()
    idx = int(index)

    if feature == 'pseudoknot':
        dotbracket = _extract_dotbracket(dssr_data)
        nts_list = dssr_data.get('nts', None)
        if nts_list is None:
            raise CmdException('No nts found in DSSR output')

        layers = parse_dotbracket_pseudoknots(dotbracket)
        if not layers:
            raise CmdException('No pseudoknot layers found')

        layer_keys = sorted(layers.keys())
        if idx < 1 or idx > len(layer_keys):
            raise CmdException('pseudoknot layer index %d out of range (1..%d)' % (idx, len(layer_keys)))

        pairs = layers[layer_keys[idx - 1]]
        sel_str = build_selection_from_layer(pairs, nts_list)
        if not sel_str:
            raise CmdException('Could not build selection for pseudoknot layer %d' % idx)
        return sel_str

    if feature not in FEATURE_MAP:
        raise CmdException('Unknown feature "%s"' % feature)

    json_key = FEATURE_MAP[feature]
    feature_list = dssr_data.get(json_key, None)
    if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
        raise CmdException('No "%s" found in DSSR output' % json_key)

    if idx < 1 or idx > len(feature_list):
        raise CmdException('%s index %d out of range (1..%d)' % (feature, idx, len(feature_list)))

    entry = feature_list[idx - 1]

    if feature == 'pairs':
        return build_selection_from_pair(entry)

    if feature == 'hbonds':
        sel_str, _, _ = build_selection_from_hbond(entry, hbonds_mode='residue')
        return sel_str

    if feature in ('stems', 'helices'):
        return build_selection_from_stem(entry)

    if feature == 'hairpins':
        return build_selection_from_hairpin(entry)

    if feature in ('stacks', 'nonstack', 'bulges', 'iloops', 'internal', 'junctions',
                   'sssegments', 'ssSegments', 'multiplets', 'splayunits'):
        nts_long = entry.get('nts_long', '')
        if not nts_long:
            raise CmdException('%s entry missing nts_long field' % feature)
        nts_list_parsed = [nt.strip() for nt in nts_long.split(',') if nt.strip()]
        return build_selection_from_nts_list(nts_list_parsed)

    if feature == 'coaxstacks':
        stems_list = dssr_data.get('stems', [])
        if not stems_list:
            raise CmdException('No stems found, required for coaxStacks')
        return build_selection_from_coaxstack(entry, stems_list)

    if feature == 'atom2bases':
        return build_selection_from_atom2base(entry)

    if feature == 'aminors':
        return build_selection_from_aminor(entry)

    if feature == 'nts':
        nt_id = entry.get('nt_id')
        if not nt_id:
            raise CmdException('Nucleotide entry missing nt_id field')
        c, r = parse_nt_id(nt_id)
        return '(chain %s and resi %s)' % (c, r)

    raise CmdException('Feature "%s" not supported for residue selection' % feature)


def dssr_select_v3(selection='all',
                   state=-1,
                   feature='pairs',
                   index=1,
                   name='dssr_select_v3',
                   exe='x3dna-dssr',
                   show_info=0,
                   quiet=1,
                   color='auto',
                   precolor=1,
                   hbonds_mode='residue',
                   distance_name=''):
    import tempfile, os

    state = int(state)
    index = int(index)
    show_info = int(show_info)
    quiet = int(quiet)
    precolor = int(precolor)
    feature = unquote(feature).lower().strip()

    user_color = _resolve_color_spec(unquote(color).strip())
    hb_mode = str(hbonds_mode).strip().lower()

    if feature in ('features', 'help'):
        keys = sorted(FEATURE_MAP.keys())
        print('Supported features: ' + ', '.join(keys))
        print('Tip: use index=0 to list detected items for a feature.')
        return

    if feature not in FEATURE_MAP:
        valid = ', '.join(sorted(FEATURE_MAP.keys()))
        raise CmdException('Unknown feature "%s". Valid: %s' % (feature, valid))

    json_key = FEATURE_MAP[feature]

    if state == 0 or state < 0:
        state = cmd.get_state()

    tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
    tmpfilepdb = tmp.name
    tmp.close()

    try:
        cmd.save(tmpfilepdb, selection, state)
        if precolor:
            cmd.color('gray', selection)

        dssr_data = run_dssr_json(tmpfilepdb, exe)

        if feature == 'pseudoknot':
            layer_colors = ['blue', 'pink', 'green', 'yellow', 'orange']

            dotbracket = _extract_dotbracket(dssr_data)
            nts_list = dssr_data.get('nts', None)
            if nts_list is None:
                raise CmdException('No nts found in DSSR output')

            layers = parse_dotbracket_pseudoknots(dotbracket)
            if not layers:
                raise CmdException('No pseudoknot layers found in structure')

            layer_keys = sorted(layers.keys())

            if index == 0:
                print('pseudoknot: %d layer(s)' % len(layer_keys))
                for j, k in enumerate(layer_keys, 1):
                    print('  layer %d (key=%s): %d pair(s)' % (j, str(k), len(layers[k])))
                if not quiet and show_info:
                    print('pseudoknot dot-bracket: ' + str(dotbracket))
                return

            if index < 1 or index > len(layer_keys):
                raise CmdException('Layer index %d out of range (1..%d)' % (index, len(layer_keys)))

            layer_key = layer_keys[index - 1]
            pairs = layers[layer_key]
            layer_color = layer_colors[(index - 1) % len(layer_colors)]

            sel_str = build_selection_from_layer(pairs, nts_list)
            if sel_str is None:
                raise CmdException('Could not build selection for layer %d' % index)

            cmd.select(name, sel_str)
            cmd.color(user_color if user_color else layer_color, name)

            if not quiet:
                print('dssr_select_v3: selection "%s" pseudoknot layer %d with %d pair(s)' % (name, index, len(pairs)))
            return

        feature_list = dssr_data.get(json_key, None)
        if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
            raise CmdException('No "%s" found in DSSR output' % json_key)

        if index == 0:
            total = len(feature_list)
            print('%s: %d item(s)' % (feature, total))
            show_n = 20 if not quiet else 10
            show_n = min(show_n, total)
            for i in range(show_n):
                print('  ' + _preview_entry(feature, feature_list[i], i + 1))
            if total > show_n:
                print('  ... (%d more)' % (total - show_n))
            return

        if index < 1 or index > len(feature_list):
            raise CmdException('%s index %d out of range (1..%d)' % (feature, index, len(feature_list)))

        entry = feature_list[index - 1]

        dist_a1 = None
        dist_a2 = None

        if feature == 'pairs':
            sel_str = build_selection_from_pair(entry)
        elif feature == 'hbonds':
            sel_str, dist_a1, dist_a2 = build_selection_from_hbond(entry, hbonds_mode=hb_mode)
        elif feature in ('stems', 'helices'):
            sel_str = build_selection_from_stem(entry)
        elif feature == 'hairpins':
            sel_str = build_selection_from_hairpin(entry)
        elif feature in ('stacks', 'nonstack', 'bulges', 'iloops', 'internal', 'junctions',
                         'sssegments', 'ssSegments', 'multiplets', 'splayunits'):
            nts_long = entry.get('nts_long', '')
            if not nts_long:
                raise CmdException('%s entry missing nts_long field' % feature)
            nts_list_parsed = [nt.strip() for nt in nts_long.split(',') if nt.strip()]
            sel_str = build_selection_from_nts_list(nts_list_parsed)
        elif feature == 'coaxstacks':
            stems_list = dssr_data.get('stems', [])
            if not stems_list:
                raise CmdException('No stems found, required for coaxStacks')
            sel_str = build_selection_from_coaxstack(entry, stems_list)
        elif feature == 'atom2bases':
            sel_str = build_selection_from_atom2base(entry)
        elif feature == 'aminors':
            sel_str = build_selection_from_aminor(entry)
        elif feature == 'nts':
            nt_id = entry.get('nt_id')
            if not nt_id:
                raise CmdException('Nucleotide entry missing nt_id field')
            c, r = parse_nt_id(nt_id)
            sel_str = '(chain %s and resi %s)' % (c, r)
        else:
            raise CmdException('Feature type "%s" not implemented' % feature)

        cmd.select(name, sel_str)
        cmd.color(user_color if user_color else 'pink', name)
        selected_features.append(feature)

        if feature == 'hbonds' and hb_mode == 'distance' and dist_a1 and dist_a2:
            dist_name = str(distance_name).strip() if distance_name else ''
            if not dist_name:
                dist_name = 'dist_%s' % name
            try:
                cmd.delete(dist_name)
            except Exception:
                pass
            cmd.distance(dist_name, dist_a1, dist_a2)
            try:
                cmd.show('dashes', dist_name)
            except Exception:
                pass
            try:
                cmd.color(user_color if user_color else 'yellow', dist_name)
            except Exception:
                pass

        if not quiet:
            print('dssr_select_v3: created selection "%s" for %s (index %d) in state %d' % (name, feature, index, state))

    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass


def _dssr_default_selection():
    objs = cmd.get_object_list('enabled')
    if len(objs) == 1:
        return objs[0]
    return 'all'


def dssr(sel=None, selection=None,
         f=None, feature='pairs',
         i=None, index=1,
         n=None, name=None,
         q=None, quiet=0,
         si=None, show_info=0,
         st=None, state=-1,
         exe='x3dna-dssr',
         color='auto',
         display=0,
         stick_radius=0.25,
         do_zoom=1,
         pc=None,
         precolor=1,
         hbmode=None,
         hbonds_mode='residue',
         distname=None,
         distance_name=''):
    selection = selection or sel or _dssr_default_selection()

    feature_in = f if f is not None else feature
    feature_in = unquote(feature_in).strip()

    if feature_in.lower() in ('features', 'help'):
        dssr_select_v3(selection=selection, state=state, feature='features', index=0,
                      name='dssr_select_v3', exe=exe, show_info=0,
                      quiet=int(q if q is not None else quiet),
                      color='auto', precolor=int(precolor))
        return

    idx = int(i if i is not None else index)
    qt = int(q if q is not None else quiet)
    si2 = int(si if si is not None else show_info)
    st2 = int(st if st is not None else state)

    if pc is not None:
        precolor = int(pc)
    precolor = int(precolor)

    hb_mode = hbmode if hbmode is not None else hbonds_mode
    hb_mode = str(hb_mode).strip().lower()

    dist_nm = distname if distname is not None else distance_name
    dist_nm = str(dist_nm).strip()

    nm = name if name is not None else n
    if not nm:
        nm = '%s%d' % (feature_in.lower(), idx)

    dssr_select_v3(selection=selection, state=st2, feature=feature_in, index=idx,
                   name=nm, exe=exe, show_info=si2, quiet=qt, color=color,
                   precolor=precolor, hbonds_mode=hb_mode, distance_name=dist_nm)

    if int(display):
        cmd.show('sticks', nm)
        try:
            cmd.set('stick_radius', float(stick_radius), nm)
        except Exception:
            pass
        if int(do_zoom):
            cmd.zoom(nm)


def _unused_name(prefix):
    try:
        return cmd.get_unused_name(prefix)
    except Exception:
        base = str(prefix) if prefix else 'obj'
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
            name = '%s%d' % (base, k)


def dssr_block(selection='all',
               state=-1,
               block_file='face',
               block_depth=0.5,
               name='',
               exe='x3dna-dssr',
               quiet=1):
    import subprocess
    import tempfile, os

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

    tmp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
    tmp_r3d = tempfile.NamedTemporaryFile(suffix='.r3d', delete=False)
    tmpfilepdb = tmp_pdb.name
    tmpfiler3d = tmp_r3d.name
    tmp_pdb.close()
    tmp_r3d.close()

    if not name:
        name = _unused_name('dssr_block')

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
                '--block-file=' + unquote(block_file),
                '--block-depth=' + str(block_depth),
                '-i=' + tmpfilepdb,
                '-o=' + tmpfiler3d,
            ]

            try:
                p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                rc = p.returncode
            except OSError:
                raise CmdException('Cannot execute exe="%s"' % exe)

            if rc != 0:
                err_txt = ''
                try:
                    err_txt = err.decode('utf-8', errors='replace') if err else ''
                except Exception:
                    err_txt = str(err)
                raise CmdException('DSSR block failed (rc=%s). stderr tail: %s' % (str(rc), _safe_tail(err_txt)))

            cmd.load(tmpfiler3d, name, max(1, st), zoom=0)

        if not quiet:
            print('dssr_block: loaded "%s" (block_file=%s, block_depth=%s)' % (name, str(block_file), str(block_depth)))

    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass
        try:
            os.remove(tmpfiler3d)
        except OSError:
            pass


def _wrap_seq(s, width):
    try:
        width = int(width)
    except Exception:
        width = 80
    if width <= 0:
        return s
    return '\n'.join(s[i:i + width] for i in range(0, len(s), width))


def _revcomp(seq):
    s = ''.join([c for c in str(seq).upper() if c.isalpha()])
    if 'U' in s and 'T' not in s:
        comp = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    else:
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in s[::-1])


def _parse_fastastr(fasta_text):
    blocks = []
    header = None
    seq = []
    for line in str(fasta_text).splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                blocks.append((header, ''.join(seq)))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        blocks.append((header, ''.join(seq)))
    return blocks


def dssr_seq(selection='all',
             chain='',
             fmt='raw',
             wrap=80,
             rc=0,
             quiet=1):
    fmt = str(fmt).strip().lower()
    quiet = int(quiet)
    rc = int(rc)

    sel = selection
    ch = str(chain).strip()
    if ch and ch.lower() != 'all':
        sel = '(%s) and chain %s' % (selection, ch)

    try:
        fasta = cmd.get_fastastr(sel)
    except Exception as e:
        raise CmdException('get_fastastr failed: %s' % e)

    blocks = _parse_fastastr(fasta)
    if not blocks:
        raise CmdException('No FASTA sequence extracted from selection="%s"' % sel)

    out_lines = []
    for hdr, seq in blocks:
        s = ''.join([c for c in seq.upper() if c.isalpha()])
        if rc:
            s = _revcomp(s)

        if fmt == 'fasta':
            out_lines.append('>' + hdr)
            out_lines.append(_wrap_seq(s, wrap))
        else:
            out_lines.append(hdr + ': ' + (_wrap_seq(s, wrap) if int(wrap) > 0 else s))

    out = '\n'.join(out_lines)
    if not quiet:
        print(out)
    return out


def _restore_pretty_colors(objs):
    for o in objs:
        try:
            cmd.spectrum('resi', 'rainbow', '(%s) and polymer' % o)
        except Exception:
            pass
        try:
            cmd.color('atomic', '(%s) and not polymer' % o)
        except Exception:
            pass


def _clear_keep_molecules(apply_gray):
    try:
        all_objs = list(cmd.get_object_list())
    except Exception:
        all_objs = []

    keep = []
    delete_objs = []
    for o in all_objs:
        try:
            n = int(cmd.count_atoms(o))
        except Exception:
            n = 0
        if n > 0:
            keep.append(o)
        else:
            delete_objs.append(o)

    for o in delete_objs:
        try:
            cmd.delete(o)
        except Exception:
            pass

    for o in list(_DSSR_BLOCK_OBJECTS):
        try:
            cmd.delete(o)
        except Exception:
            pass
    _DSSR_BLOCK_OBJECTS.clear()

    try:
        for m in cmd.get_names('measurements'):
            try:
                cmd.delete(m)
            except Exception:
                pass
    except Exception:
        pass

    try:
        for s in cmd.get_names('selections'):
            try:
                cmd.delete(s)
            except Exception:
                pass
    except Exception:
        pass

    try:
        cmd.label('all', '')
    except Exception:
        pass

    try:
        cmd.select('sele', 'none')
    except Exception:
        pass

    if apply_gray:
        for o in keep:
            try:
                cmd.color('gray', o)
            except Exception:
                pass
    else:
        _restore_pretty_colors(keep)

    try:
        cmd.zoom('all')
    except Exception:
        pass
    try:
        cmd.reset()
    except Exception:
        pass


try:
    from pymol.Qt import QtWidgets, QtCore
except Exception:
    QtWidgets = None
    QtCore = None

_DSSR_GUI_DIALOG = None

class _DSSRGuiDialog(QtWidgets.QDialog if QtWidgets else object):
    def __init__(self):
        super(_DSSRGuiDialog, self).__init__()
        self.setWindowTitle('DSSR GUI')
        self.resize(1020, 690)

        self._current_feature = 'pairs'
        self._cache_key = None
        self._cache_data = None
        self._items_all = []
        self._items_filtered = []
        self._page = 0

        root = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QGridLayout()
        root.addLayout(top)

        self.obj_combo = QtWidgets.QComboBox()
        self.obj_combo.currentIndexChanged.connect(self._on_object_changed)

        self.state_combo = QtWidgets.QComboBox()
        self.count_btn = QtWidgets.QPushButton('count')
        self.count_btn.clicked.connect(self._count_states_clicked)

        self.exe_edit = QtWidgets.QLineEdit('x3dna-dssr')

        self.color_edit = QtWidgets.QLineEdit('auto')
        self.name_edit = QtWidgets.QLineEdit('')

        self.hb_mode_combo = QtWidgets.QComboBox()
        self.hb_mode_combo.addItems(['residue', 'atom', 'distance'])

        self.precolor_cb = QtWidgets.QCheckBox('gray precolor')
        self.precolor_cb.setChecked(True)

        self.display_cb = QtWidgets.QCheckBox('display sticks')
        self.display_cb.setChecked(False)

        self.zoom_cb = QtWidgets.QCheckBox('zoom')
        self.zoom_cb.setChecked(True)

        self.showinfo_cb = QtWidgets.QCheckBox('show_info (pseudoknot)')
        self.showinfo_cb.setChecked(False)

        self.radius_spin = QtWidgets.QDoubleSpinBox()
        self.radius_spin.setMinimum(0.01)
        self.radius_spin.setMaximum(5.0)
        self.radius_spin.setSingleStep(0.05)
        self.radius_spin.setValue(0.25)

        self.rna_btn = QtWidgets.QPushButton('RNA only')
        self.rna_btn.clicked.connect(self._make_rna_only)

        self.status_label = QtWidgets.QLabel('')
        self.status_label.setWordWrap(True)

        self.block_file_combo = QtWidgets.QComboBox()
        self.block_file_combo.setEditable(True)
        self.block_file_combo.addItems(['face', 'edge', 'wc', 'equal', 'minor', 'gray', 'wc-minor'])

        self.block_depth_spin = QtWidgets.QDoubleSpinBox()
        self.block_depth_spin.setMinimum(0.01)
        self.block_depth_spin.setMaximum(5.0)
        self.block_depth_spin.setSingleStep(0.05)
        self.block_depth_spin.setValue(0.5)

        self.make_blocks_btn = QtWidgets.QPushButton('make blocks')
        self.make_blocks_btn.clicked.connect(self._make_blocks_clicked)

        self.seq_btn = QtWidgets.QPushButton('seq view')
        self.seq_btn.clicked.connect(self._seq_view_clicked)

        top.addWidget(QtWidgets.QLabel('object'), 0, 0)
        top.addWidget(self.obj_combo, 0, 1, 1, 2)
        top.addWidget(self.rna_btn, 0, 3)

        top.addWidget(QtWidgets.QLabel('state'), 0, 4)
        top.addWidget(self.state_combo, 0, 5)
        top.addWidget(self.count_btn, 0, 6)

        top.addWidget(QtWidgets.QLabel('exe'), 1, 0)
        top.addWidget(self.exe_edit, 1, 1, 1, 6)

        top.addWidget(QtWidgets.QLabel('color'), 2, 0)
        top.addWidget(self.color_edit, 2, 1)
        top.addWidget(QtWidgets.QLabel('name'), 2, 2)
        top.addWidget(self.name_edit, 2, 3)
        top.addWidget(QtWidgets.QLabel('hbonds_mode'), 2, 4)
        top.addWidget(self.hb_mode_combo, 2, 5, 1, 2)

        top.addWidget(QtWidgets.QLabel('block_file'), 3, 0)
        top.addWidget(self.block_file_combo, 3, 1)
        top.addWidget(QtWidgets.QLabel('block_depth'), 3, 2)
        top.addWidget(self.block_depth_spin, 3, 3)
        top.addWidget(self.make_blocks_btn, 3, 4, 1, 3)

        opts = QtWidgets.QHBoxLayout()
        opts.addWidget(self.precolor_cb)
        opts.addWidget(self.display_cb)
        opts.addWidget(self.zoom_cb)
        opts.addWidget(self.showinfo_cb)
        opts.addWidget(QtWidgets.QLabel('stick_radius'))
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
            b = QtWidgets.QPushButton(feat)
            b.setCheckable(True)
            b.clicked.connect(self._make_feature_handler(feat))
            r = idx // cols
            c = idx % cols
            btn_grid.addWidget(b, r, c)
            self._feature_buttons[feat] = b
        if 'pairs' in self._feature_buttons:
            self._feature_buttons['pairs'].setChecked(True)

        root.addWidget(btn_area)

        bar = QtWidgets.QHBoxLayout()

        self.filter_edit = QtWidgets.QLineEdit('')
        self.filter_edit.setPlaceholderText('filter...')
        self.filter_edit.textChanged.connect(self._on_filter_changed)

        self.page_size_spin = QtWidgets.QSpinBox()
        self.page_size_spin.setMinimum(50)
        self.page_size_spin.setMaximum(5000)
        self.page_size_spin.setValue(500)
        self.page_size_spin.valueChanged.connect(self._on_page_size_changed)

        self.prev_btn = QtWidgets.QPushButton('Prev')
        self.next_btn = QtWidgets.QPushButton('Next')
        self.prev_btn.clicked.connect(self._prev_page)
        self.next_btn.clicked.connect(self._next_page)

        self.page_label = QtWidgets.QLabel('')

        self.refresh_obj_btn = QtWidgets.QPushButton('refresh objects')
        self.refresh_obj_btn.clicked.connect(self.refresh_objects)
        self.refresh_list_btn = QtWidgets.QPushButton('refresh list')
        self.refresh_list_btn.clicked.connect(self.refresh_list)

        self.zoom_all_btn = QtWidgets.QPushButton('zoom all')
        self.zoom_all_btn.clicked.connect(self._zoom_all)

        self.reset_view_btn = QtWidgets.QPushButton('reset view')
        self.reset_view_btn.clicked.connect(self._reset_view)

        bar.addWidget(QtWidgets.QLabel('filter'))
        bar.addWidget(self.filter_edit, 2)
        bar.addWidget(QtWidgets.QLabel('max/page'))
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

        report_bar = QtWidgets.QHBoxLayout()
        self.report_btn = QtWidgets.QPushButton('Generate RNA Report')
        self.report_btn.clicked.connect(self._generate_rna_report_clicked)

        self.make_all_sel_btn = QtWidgets.QPushButton('make selections (all)')
        self.make_all_sel_btn.clicked.connect(self._make_all_selections_clicked)

        self.auto_color_btn = QtWidgets.QPushButton('auto color')
        self.auto_color_btn.clicked.connect(self._auto_color_clicked)

        self.clear_report_btn = QtWidgets.QPushButton('clear report')
        self.clear_report_btn.clicked.connect(self._clear_report)

        report_bar.addWidget(self.report_btn)
        report_bar.addWidget(self.make_all_sel_btn)
        report_bar.addWidget(self.auto_color_btn)
        report_bar.addWidget(self.clear_report_btn)
        report_bar.addStretch(1)
        root.addLayout(report_bar)
        self.list_widget = QtWidgets.QListWidget()
        try:
            self.list_widget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        except Exception:
            pass
        self.list_widget.itemClicked.connect(self._on_item_clicked_preview)
        self.list_widget.itemDoubleClicked.connect(self._on_item_double_clicked)
        self.report_box = QtWidgets.QPlainTextEdit()
        self.report_box.setReadOnly(True)
        try:
            self.report_box.setPlaceholderText('RNA Structure Summary will appear here...')
        except Exception:
            pass

        split = QtWidgets.QSplitter()
        split.setOrientation(QtCore.Qt.Horizontal)
        split.addWidget(self.list_widget)
        split.addWidget(self.report_box)
        try:
            split.setStretchFactor(0, 3)
            split.setStretchFactor(1, 2)
        except Exception:
            pass

        root.addWidget(split, 1)

        self.refresh_objects()
        self._update_state_combo()
        self.refresh_list()

    def _enable_seq_view(self):
        try:
            cmd.set('seq_view', 1)
        except Exception:
            try:
                cmd.do('set seq_view, 1')
            except Exception:
                pass

    def _seq_view_clicked(self):
        self._enable_seq_view()
        sel = self._get_object_text()
        try:
            cmd.select('sele', '(%s) and polymer.nucleic' % sel)
        except Exception:
            pass

    def _make_feature_handler(self, feat):
        def handler():
            self._current_feature = feat
            for k, b in self._feature_buttons.items():
                b.setChecked(k == feat)
            self.refresh_list()
        return handler

    def _on_object_changed(self):
        self._update_state_combo()

    def _count_states_clicked(self):
        obj = self._get_object_text()
        try:
            n = int(cmd.count_states(obj))
        except Exception:
            n = 0
        self._update_state_combo(force_n=n)
        try:
            QtWidgets.QMessageBox.information(self, 'count_states', 'count_states %s = %d' % (obj, n))
        except Exception:
            pass

    def _get_object_text(self):
        txt = self.obj_combo.currentText().strip() if self.obj_combo.count() else ''
        return txt if txt else 'all'

    def _update_state_combo(self, force_n=None):
        obj = self._get_object_text()
        try:
            n = int(force_n) if force_n is not None else int(cmd.count_states(obj))
        except Exception:
            n = 0

        current = self.state_combo.currentData() if self.state_combo.count() else None
        self.state_combo.clear()

        self.state_combo.addItem('current', -1)
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
        try:
            objs = cmd.get_object_list('enabled')
            if not objs:
                objs = cmd.get_object_list()
        except Exception:
            objs = []

        current = self.obj_combo.currentText().strip() if self.obj_combo.count() else ''
        self.obj_combo.clear()

        if not objs:
            self.obj_combo.addItem('all')
        else:
            for o in objs:
                self.obj_combo.addItem(o)

        if current:
            i = self.obj_combo.findText(current)
            if i >= 0:
                self.obj_combo.setCurrentIndex(i)
        else:
            try:
                enabled = cmd.get_object_list('enabled')
            except Exception:
                enabled = []
            if len(enabled) == 1:
                i = self.obj_combo.findText(enabled[0])
                if i >= 0:
                    self.obj_combo.setCurrentIndex(i)

        self._update_state_combo()

    def _zoom_all(self):
        try:
            cmd.zoom('all')
        except Exception:
            pass

    def _reset_view(self):
        apply_gray = True if self.precolor_cb.isChecked() else False
        _clear_keep_molecules(apply_gray)
        self.refresh_objects()
        self.refresh_list()

    def _unique_object_name(self, base):
        base = str(base)
        if not base:
            base = 'rna_only'
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
            name = '%s%d' % (base, k)

    def _make_rna_only(self):
        obj = self._get_object_text()
        new_name = self._unique_object_name('rna_only')
        try:
            cmd.create(new_name, '(%s) and polymer.nucleic' % obj)
        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, 'RNA only', str(e))
            except Exception:
                pass
            return

        self.refresh_objects()
        i = self.obj_combo.findText(new_name)
        if i >= 0:
            self.obj_combo.setCurrentIndex(i)
        self.refresh_list()

    def _clear_report(self):
        try:
            self.report_box.setPlainText('')
        except Exception:
            pass

    def _get_dssr_context(self):
        sel_obj = self._get_object_text()
        exe = self.exe_edit.text().strip() or 'x3dna-dssr'
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0
        return sel_obj, exe, st, precolor_on

    def _append_report(self, text):
        try:
            cur = self.report_box.toPlainText()
        except Exception:
            cur = ''
        if cur:
            out = cur.rstrip('\n') + '\n\n' + str(text).rstrip('\n') + '\n'
        else:
            out = str(text).rstrip('\n') + '\n'
        try:
            self.report_box.setPlainText(out)
        except Exception:
            pass

    def _generate_rna_report_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            report = _format_rna_summary_text(dssr_data)
            try:
                self.report_box.setPlainText(report + '\n')
            except Exception:
                self._append_report(report)
        except Exception as e:
            msg = 'RNA report error: %s' % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, 'RNA Report', msg)
            except Exception:
                pass

    def _make_all_selections_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            made, skipped = self._make_all_selections(dssr_data, sel_obj)
            self._append_report('Created selections: %s' % (' '.join(made) if made else '(none)'))
            if skipped:
                self._append_report('Skipped (empty): %s' % (' '.join(skipped)))
        except Exception as e:
            msg = 'make selections error: %s' % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, 'make selections', msg)
            except Exception:
                pass

    def _auto_color_clicked(self):
        sel_obj, exe, st, precolor_on = self._get_dssr_context()
        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            self._make_all_selections(dssr_data, sel_obj)
            self._apply_auto_colors()
            self._append_report('Auto color applied: stems green, hairpins blue, pseudoknots red, aminors purple')
        except Exception as e:
            msg = 'auto color error: %s' % str(e)
            self._append_report(msg)
            try:
                QtWidgets.QMessageBox.critical(self, 'auto color', msg)
            except Exception:
                pass

    def _make_all_selections(self, dssr_data, obj_sel):
        targets = [
            ('pairs', 'pairs_all'),
            ('hairpins', 'hairpins_all'),
            ('stems', 'stems_all'),
            ('bulges', 'bulges_all'),
            ('junctions', 'junctions_all'),
            ('pseudoknot', 'pseudoknots_all'),
            ('aminors', 'aminors_all'),
            ('stacks', 'stacks_all'),
        ]

        made = []
        skipped = []

        for feat, name in targets:
            residues = _collect_residues_all(dssr_data, feat)
            sel_core = _compact_sel_from_residues(residues)
            if not sel_core:
                skipped.append(name)
                try:
                    cmd.delete(name)
                except Exception:
                    pass
                continue

            expr = '((%s) and polymer.nucleic and (%s))' % (obj_sel, sel_core)
            try:
                cmd.select(name, expr)
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
            cmd.color('green', 'stems_all')
        except Exception:
            pass
        try:
            cmd.color('blue', 'hairpins_all')
        except Exception:
            pass
        try:
            cmd.color('red', 'pseudoknots_all')
        except Exception:
            pass
        try:
            cmd.color('purple', 'aminors_all')
        except Exception:
            pass


    def _big_object_warning(self, sel):
        thresh = 250000
        try:
            n_atoms = int(cmd.count_atoms(sel))
        except Exception:
            n_atoms = 0
        if n_atoms >= thresh:
            self.status_label.setText('warning: large selection (%d atoms). recommended: use RNA only and avoid feature=nts/hbonds first.' % n_atoms)
        else:
            self.status_label.setText('')

    def _get_dssr_data(self, selection, state, exe, precolor_on):
        import tempfile, os
        cache_key = (str(selection), int(state), str(exe))
        if self._cache_key == cache_key and self._cache_data is not None:
            if int(precolor_on):
                try:
                    cmd.color('gray', selection)
                except Exception:
                    pass
            return self._cache_data

        tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
        tmpfilepdb = tmp.name
        tmp.close()

        try:
            cmd.save(tmpfilepdb, selection, state)
            if int(precolor_on):
                cmd.color('gray', selection)
            data = run_dssr_json(tmpfilepdb, exe)
        finally:
            try:
                os.remove(tmpfilepdb)
            except OSError:
                pass

        self._cache_key = cache_key
        self._cache_data = data
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
        exe = self.exe_edit.text().strip() or 'x3dna-dssr'
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0

        self._big_object_warning(sel)

        self.list_widget.clear()
        self.list_widget.addItem('loading...')
        QtWidgets.QApplication.processEvents()

        try:
            dssr_data = self._get_dssr_data(sel, st, exe, precolor_on)

            items = []

            if feat == 'pseudoknot':
                dotbracket = _extract_dotbracket(dssr_data)
                nts_list = dssr_data.get('nts', None)
                if nts_list is None:
                    raise CmdException('No nts found in DSSR output')

                layers = parse_dotbracket_pseudoknots(dotbracket)
                if not layers:
                    raise CmdException('No pseudoknot layers found')

                layer_keys = sorted(layers.keys())
                for j, k in enumerate(layer_keys, 1):
                    line = '%d: layer key=%s pairs=%d' % (j, str(k), len(layers[k]))
                    items.append((j, line))

            else:
                if feat not in FEATURE_MAP:
                    raise CmdException('Unknown feature "%s"' % feat)
                json_key = FEATURE_MAP[feat]
                feature_list = dssr_data.get(json_key, None)
                if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
                    raise CmdException('No "%s" found in DSSR output' % json_key)

                total = len(feature_list)
                for i in range(total):
                    line = _preview_entry(feat, feature_list[i], i + 1)
                    items.append((i + 1, line))

            self._items_all = items
            self._page = 0
            self._render_list()

        except Exception as e:
            self._items_all = []
            self._items_filtered = []
            self.list_widget.clear()
            self.list_widget.addItem('ERROR: %s' % str(e))
            try:
                print('dssr_gui error: %s' % str(e))
            except Exception:
                pass

    def _render_list(self):
        q = self.filter_edit.text().strip().lower()
        if q:
            self._items_filtered = [it for it in self._items_all if q in it[1].lower()]
        else:
            self._items_filtered = list(self._items_all)

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

        self.page_label.setText('items %d, filtered %d, page %d/%d' % (total_all, total_f, self._page + 1, pages))
        self.prev_btn.setEnabled(self._page > 0)
        self.next_btn.setEnabled(self._page + 1 < pages)

        try:
            self.setWindowTitle('DSSR GUI - %s' % self._current_feature)
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
        exe = self.exe_edit.text().strip() or 'x3dna-dssr'
        st = self._get_state_value()
        precolor_on = 1 if self.precolor_cb.isChecked() else 0

        try:
            dssr_data = self._get_dssr_data(sel_obj, st, exe, precolor_on)
            sel_str = _build_residue_sel_from_dssr(dssr_data, feat, idx)
            cmd.select('sele', sel_str)
            cmd.select('sele', 'byres (sele) and polymer.nucleic')
        except Exception as e:
            try:
                self.status_label.setText('preview error: %s' % str(e))
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
        exe = self.exe_edit.text().strip() or 'x3dna-dssr'
        st = self._get_state_value()

        nm = self.name_edit.text().strip()
        if not nm:
            nm = '%s%d' % (feat, idx)

        col = self.color_edit.text().strip() or 'auto'
        precolor_on = 1 if self.precolor_cb.isChecked() else 0
        display_on = 1 if self.display_cb.isChecked() else 0
        zoom_on = 1 if self.zoom_cb.isChecked() else 0
        showinfo_on = 1 if self.showinfo_cb.isChecked() else 0
        radius = float(self.radius_spin.value())

        hb_mode = self.hb_mode_combo.currentText().strip().lower()
        dist_name = 'dist_%s' % nm

        try:
            dssr(sel=sel, f=feat, i=idx, n=nm, q=0, si=showinfo_on, st=st,
                 exe=exe, color=col, display=display_on, stick_radius=radius,
                 do_zoom=zoom_on, pc=precolor_on, hbmode=hb_mode, distname=dist_name)
        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, 'DSSR GUI error', str(e))
            except Exception:
                pass
            try:
                print('dssr_gui select error: %s' % str(e))
            except Exception:
                pass

    def _make_blocks_clicked(self):
        self._enable_seq_view()

        sel_obj = self._get_object_text()
        exe = self.exe_edit.text().strip() or 'x3dna-dssr'
        st = self._get_state_value()

        block_file = self.block_file_combo.currentText().strip() or 'face'
        block_depth = float(self.block_depth_spin.value())

        base = self.name_edit.text().strip()
        if not base:
            base = 'blk'

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
                    QtWidgets.QMessageBox.critical(self, 'make blocks error', str(e))
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

                    sel_str = _build_residue_sel_from_dssr(dssr_data, feat, idx)
                    sel_for_block = 'byres (%s) and polymer.nucleic' % sel_str

                    obj_name = '%s_%s_%d' % (base, feat, idx)
                    try:
                        cmd.delete(obj_name)
                    except Exception:
                        pass

                    dssr_block(selection=sel_for_block, state=st,
                               block_file=block_file, block_depth=block_depth,
                               name=obj_name, exe=exe, quiet=1)

                    _DSSR_BLOCK_OBJECTS.add(obj_name)
                    made_names.append(obj_name)

            else:
                try:
                    n = int(cmd.count_atoms('sele'))
                except Exception:
                    n = 0
                if n > 0:
                    sel_for_block = 'byres (sele) and polymer.nucleic'
                else:
                    sel_for_block = '(%s) and polymer.nucleic' % sel_obj

                obj_name = '%s_sel' % base
                try:
                    cmd.delete(obj_name)
                except Exception:
                    pass

                dssr_block(selection=sel_for_block, state=st,
                           block_file=block_file, block_depth=block_depth,
                           name=obj_name, exe=exe, quiet=1)

                _DSSR_BLOCK_OBJECTS.add(obj_name)
                made_names.append(obj_name)

            if made_names and self.zoom_cb.isChecked():
                try:
                    cmd.zoom('(%s)' % ' or '.join(made_names))
                except Exception:
                    try:
                        cmd.zoom('all')
                    except Exception:
                        pass

        except Exception as e:
            try:
                QtWidgets.QMessageBox.critical(self, 'make blocks error', str(e))
            except Exception:
                pass
            try:
                print('make blocks error: %s' % str(e))
            except Exception:
                pass


def dssr_gui():
    global _DSSR_GUI_DIALOG
    if QtWidgets is None or QtCore is None:
        raise CmdException('Qt is not available in this PyMOL build')
    if _DSSR_GUI_DIALOG is None:
        _DSSR_GUI_DIALOG = _DSSRGuiDialog()
    _DSSR_GUI_DIALOG.show()
    _DSSR_GUI_DIALOG.raise_()
    _DSSR_GUI_DIALOG.activateWindow()


cmd.extend('dssr_select_v3', dssr_select_v3)
cmd.extend('dssr_gui', dssr_gui)
cmd.extend('dssr_block', dssr_block)
cmd.extend('dssr_seq', dssr_seq)

cmd.auto_arg[0].update({'dssr_select_v3': cmd.auto_arg[0]['zoom']})
cmd.auto_arg[0].update({'dssr_gui': cmd.auto_arg[0]['zoom']})
cmd.auto_arg[0].update({'dssr_block': cmd.auto_arg[0]['zoom']})
cmd.auto_arg[0].update({'dssr_seq': cmd.auto_arg[0]['zoom']})

try:
    print('Loaded DSSR helper %s' % __DSSR_PLUGIN_VERSION__)
except Exception:
    pass
