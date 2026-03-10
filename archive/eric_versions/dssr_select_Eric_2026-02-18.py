# DSSR integration with PyMOL (selection helper)
# Updated: 2026-02-18
#run /home/eric/dssr_select_v3copy.py


from pymol import cmd, CmdException


def unquote(s):
    s = str(s)
    if s.rstrip()[-1:] not in ('"', "'"):
        return s
    return cmd.safe_eval(s)


def run_dssr_json(pdb_path, exe):
    import subprocess, json
    args = [exe, '--json', '-i=' + pdb_path]
    try:
        out = subprocess.check_output(args)
    except subprocess.CalledProcessError:
        raise CmdException('"%s" failed' % exe)
    except OSError:
        raise CmdException('Cannot execute exe="%s"' % exe)

    try:
        text = out.decode('utf-8')
        data = json.loads(text)
    except Exception as e:
        raise CmdException('Failed to parse DSSR JSON: %s' % e)
    return data


def parse_nt_id(nt_id):
    """
    Robust DSSR nt id parsing.
    Accepts:
      A.U7
      A.5MC49
      A.H2U16
      ///A/145
    Returns (chain, resi) where resi is the trailing integer (may include '-').
    """
    import re
    s = str(nt_id).strip()

    # dot format: A.H2U16
    if '.' in s and not s.startswith('///'):
        try:
            chain, res = s.split('.', 1)
        except ValueError:
            raise CmdException('Unexpected nt_id format: "%s"' % s)

        m = re.search(r'(-?\d+)$', res)
        if not m:
            raise CmdException('Could not extract residue number from "%s"' % res)
        return chain, m.group(1)

    # slash format: ///A/145
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
    # used for hbonds: something@A.U7 like patterns are handled by split('@',1)
    try:
        _, rest = atom_id.split('@', 1)
    except ValueError:
        raise CmdException('Unexpected atom_id format: "%s"' % atom_id)
    return parse_nt_id(rest)


def parse_a2b_atom(atom_id):
    """
    atom2bases atom formats seen in DSSR JSON:
      OP2@A.H2U16
      O4'@A.U59
      ///A/147/O4'
    Returns (chain, resi, atom_name)
    """
    import re
    s = str(atom_id).strip()

    # format: OP2@A.H2U16
    if '@' in s:
        atom_name, nt = s.split('@', 1)
        chain, resi = parse_nt_id(nt)
        return chain, resi, atom_name

    # format: ///A/147/O4'
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


def build_selection_from_hbond(hb_entry):
    atom1_id = hb_entry.get('atom1_id')
    atom2_id = hb_entry.get('atom2_id')
    if not atom1_id or not atom2_id:
        raise CmdException('H-bond entry missing atom1_id or atom2_id')
    residues = {parse_atom_id(atom1_id), parse_atom_id(atom2_id)}
    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)


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
        atom_name = atom_name.replace('"', '\\"')
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

    if feature in ('stacks', 'nonstack', 'hairpins', 'bulges', 'iloops', 'internal', 'junctions', 'sssegments', 'ssSegments', 'multiplets', 'splayunits'):
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


selected_features = []


def dssr_select_v3(selection='all',
                   state=-1,
                   feature='pairs',
                   index=1,
                   name='dssr_select_v3',
                   exe='x3dna-dssr',
                   show_info=0,
                   quiet=1):

    import tempfile, os

    state = int(state)
    index = int(index)
    show_info = int(show_info)
    quiet = int(quiet)
    feature = unquote(feature).lower()

    feature_map = {
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

    if feature not in feature_map:
        valid = ', '.join(sorted(feature_map.keys()))
        raise CmdException('Unknown feature "%s". Valid: %s' % (feature, valid))

    json_key = feature_map[feature]

    if state == 0 or state < 0:
        state = cmd.get_state()

    tmpfilepdb = tempfile.mktemp('.pdb')

    try:
        cmd.save(tmpfilepdb, selection, state)
        cmd.color('gray', selection)
        dssr_data = run_dssr_json(tmpfilepdb, exe)

        if feature == 'pseudoknot':
            layer_colors = ['blue', 'pink', 'green', 'yellow', 'orange']

            if 'dbn' not in dssr_data:
                raise CmdException('No dot-bracket notation found in DSSR output')

            dbn_data = dssr_data['dbn']

            if isinstance(dbn_data, dict):
                if 'all_chains' in dbn_data and 'sstr' in dbn_data['all_chains']:
                    dotbracket = dbn_data['all_chains']['sstr']
                elif 'sstr' in dbn_data:
                    dotbracket = dbn_data['sstr']
                else:
                    dotbracket = None
                    for _, v in dbn_data.items():
                        if isinstance(v, dict) and 'sstr' in v:
                            dotbracket = v['sstr']
                            break
                    if dotbracket is None:
                        raise CmdException('Could not find sstr field in dbn data')
            else:
                dotbracket = dbn_data

            if 'nts' not in dssr_data:
                raise CmdException('No nucleotide information found in DSSR output')

            nts_list = dssr_data['nts']
            layers = parse_dotbracket_pseudoknots(dotbracket)

            if not layers:
                raise CmdException('No pseudoknot layers found in structure')

            if index not in layers:
                available = ', '.join(str(i) for i in sorted(layers.keys()))
                raise CmdException('Layer index %d not found. Available layers: %s' % (index, available))

            pairs = layers[index]
            color = layer_colors[index % len(layer_colors)]
            sel_str = build_selection_from_layer(pairs, nts_list)
            if sel_str is None:
                raise CmdException('Could not build selection for layer %d' % index)

            cmd.select(name, sel_str, state=state)
            cmd.color(color, name)
            if not quiet:
                print('dssr_select_v3: selection "%s" pseudoknot layer %d with %d pair(s)' % (name, index, len(pairs)))
            return

        feature_list = dssr_data.get(json_key, None)
        if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
            raise CmdException('No "%s" found in DSSR output' % json_key)

        # list mode
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

        if feature == 'pairs':
            sel_str = build_selection_from_pair(entry)
        elif feature == 'hbonds':
            sel_str = build_selection_from_hbond(entry)
        elif feature in ('stems', 'helices'):
            sel_str = build_selection_from_stem(entry)
        elif feature == 'hairpins':
            sel_str = build_selection_from_hairpin(entry)
        elif feature in ('stacks', 'nonstack', 'bulges', 'iloops', 'internal', 'junctions', 'sssegments', 'ssSegments', 'multiplets', 'splayunits'):
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
            raise CmdException('Feature type "%s" not yet implemented' % feature)

        cmd.select(name, sel_str, state=state)
        cmd.color('pink', name)
        selected_features.append(feature)

        if not quiet:
            print('dssr_select_v3: created selection "%s" for %s (index %d) in state %d' % (name, feature, index, state))

    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass


cmd.extend('dssr_select_v3', dssr_select_v3)
cmd.auto_arg[0].update({'dssr_select_v3': cmd.auto_arg[0]['zoom']})
