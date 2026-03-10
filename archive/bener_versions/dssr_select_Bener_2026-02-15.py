# Bener Dulger
# Date Updated: 02/15/26

'''
This is the shortened version of dssr_select, while also the first edition to be 
officially published to github. It is a plugin block to be used with PyMol to allow
for selection of visual structures of RNA/DNA. Its input is the dssr executable, such that
it then uses parsing to obtain the desired features.
'''

from pymol import cmd, CmdException
def unquote(s):
    s = str(s)
    return s if s.rstrip()[-1:] not in ('"', "'") else cmd.safe_eval(s)

def run_dssr_json(pdb_path, exe):
    import subprocess, json
    try:
        out = subprocess.check_output([exe, '--json', '-i=' + pdb_path])
        return json.loads(out.decode('utf-8'))
    except subprocess.CalledProcessError:
        raise CmdException('"%s" failed' % exe)
    except OSError:
        raise CmdException('Cannot execute exe="%s"' % exe)
    except Exception as e:
        raise CmdException('Failed to parse DSSR JSON: %s' % e)

def parse_nt_id(nt_id):
    try:
        chain, res = nt_id.split('.', 1)
        resi = ''.join(ch for ch in res if ch.isdigit() or ch == '-')
        if not resi:
            raise ValueError
        return chain, resi
    except ValueError:
        raise CmdException('Unexpected nt_id format: "%s"' % nt_id)

def parse_atom_id(atom_id):
    try:
        return parse_nt_id(atom_id.split('@', 1)[1])
    except (ValueError, IndexError):
        raise CmdException('Unexpected atom_id format: "%s"' % atom_id)

def parse_dotbracket_pseudoknots(dotbracket):
    stack_map = {'(': ')', '[': ']', '{': '}', '<': '>'}
    close_to_open = {v: k for k, v in stack_map.items()}
    stacks, layers = {}, {}
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
            if bracket_type in stacks and stacks[bracket_type]:
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

def build_selection(residues):
    return ' or '.join('(chain %s and resi %s)' % (c, r) for c, r in residues)

def build_selection_from_layer(layer_pairs, nts_list):
    residues = set()
    for open_idx, close_idx in layer_pairs:
        for idx in (open_idx, close_idx):
            if idx < len(nts_list) and nts_list[idx].get('nt_id'):
                residues.add(parse_nt_id(nts_list[idx]['nt_id']))
    return build_selection(residues) if residues else None

def build_selection_from_pair(pair_entry):
    return build_selection({parse_nt_id(pair_entry[k]) for k in ('nt1', 'nt2') if pair_entry.get(k)})

def build_selection_from_hbond(hb_entry):
    return build_selection({parse_atom_id(hb_entry[k]) for k in ('atom1_id', 'atom2_id') if hb_entry.get(k)})

def build_selection_from_nts_list(nts_list):
    if not nts_list:
        raise CmdException('Empty nucleotide list')
    return build_selection({parse_nt_id(nt_id) for nt_id in nts_list})

def build_selection_from_stem(stem_entry):
    residues = set()
    for pair in stem_entry.get('pairs', []):
        for key in ('nt1', 'nt2'):
            if pair.get(key):
                residues.add(parse_nt_id(pair[key]))
    if not residues:
        raise CmdException('Stem has no pairs')
    return build_selection(residues)

def build_selection_from_hairpin(hairpin_entry):
    nts_long = hairpin_entry.get('nts_long')
    if not nts_long:
        raise CmdException('Hairpin missing nts_long field')
    return build_selection_from_nts_list([nt.strip() for nt in nts_long.split(',')])

selected_features = []

def dssr_select(selection='all', state=-1, feature='pairs', index=1, 
                   name='dssr_select', exe='x3dna-dssr', show_info=0, quiet=1):
    import tempfile, os
    state, index, show_info, quiet = int(state), int(index), int(show_info), int(quiet)
    feature = unquote(feature).lower()
    
    feature_map = {'pairs': 'pairs', 'hbonds': 'hbonds', 'stems': 'stems', 'helices': 'helices',
                   'hairpins': 'hairpins', 'bulges': 'bulges', 'iloops': 'iloops', 'internal': 'iloops',
                   'junctions': 'junctions', 'ssSegments': 'ssSegments', 'sssegments': 'ssSegments',
                   'multiplets': 'multiplets', 'nts': 'nts', 'pseudoknot': 'pseudoknot'}
    
    if feature not in feature_map:
        raise CmdException('Unknown feature "%s". Valid: %s' % (feature, ', '.join(sorted(feature_map.keys()))))
    if state <= 0:
        state = cmd.get_state()
    tmpfilepdb = tempfile.mktemp('.pdb')
    try:
        cmd.save(tmpfilepdb, selection, state)
        cmd.color('gray', selection)
        dssr_data = run_dssr_json(tmpfilepdb, exe)
        if not quiet:
            print("DSSR top-level JSON keys:", sorted(dssr_data.keys()))
        if feature == 'pseudoknot':
            layer_colors = ['blue', 'pink', 'green', 'yellow', 'orange']
            if 'dbn' not in dssr_data:
                raise CmdException('No dot-bracket notation found in DSSR output')
            dbn_data = dssr_data['dbn']
            if isinstance(dbn_data, dict):
                dotbracket = (dbn_data.get('all_chains', {}).get('sstr') or 
                             dbn_data.get('sstr') or 
                             next((v['sstr'] for v in dbn_data.values() if isinstance(v, dict) and 'sstr' in v), None))
                if not dotbracket:
                    raise CmdException('Could not find sstr field in dbn data')
            else:
                dotbracket = dbn_data
            if not quiet:
                print('Dot-bracket notation:', dotbracket)
            if 'nts' not in dssr_data:
                raise CmdException('No nucleotide information found in DSSR output')
            layers = parse_dotbracket_pseudoknots(dotbracket)
            
            if not quiet:
                print('Layers found:', sorted(layers.keys()))
                for layer_id in sorted(layers.keys()):
                    print('  Layer %d: %d pairs' % (layer_id, len(layers[layer_id])))
            if not layers or index not in layers:
                available = ', '.join(str(i) for i in sorted(layers.keys()))
                raise CmdException('Layer index %d not found. Available layers: %s' % (index, available))
            pairs = layers[index]
            color = layer_colors[index % len(layer_colors)]
            
            if show_info:
                layer_type = 'Standard Watson-Crick pairs' if index == 0 else 'Pseudoknot layer %d' % index
                print('\n=== Pseudoknot Layer %d ===' % index)
                print('Dot-bracket notation: %s' % dotbracket)
                print('\n%s (%s): %d base pairs' % (layer_type, color, len(pairs)))
                print('\nBase pairs:')
                for open_idx, close_idx in pairs:
                    if open_idx < len(dssr_data['nts']) and close_idx < len(dssr_data['nts']):
                        print('  %s - %s' % (dssr_data['nts'][open_idx].get('nt_id', 'unknown'),
                                             dssr_data['nts'][close_idx].get('nt_id', 'unknown')))
                print('\n' + '='*50 + '\n')
            
            sel_str = build_selection_from_layer(pairs, dssr_data['nts'])
            if sel_str is None:
                raise CmdException('Could not build selection for layer %d' % index)
            cmd.select(name, sel_str, state=state)
            cmd.color(color, name)
            if not quiet:
                layer_desc = 'standard pairs' if index == 0 else 'pseudoknot layer %d' % index
                print('dssr_select_v3: created selection "%s" for %s (index %d) with %d base pairs (%s) in state %d' 
                      % (name, layer_desc, index, len(pairs), color, state))
        else:
            json_key = feature_map[feature]
            feature_list = dssr_data.get(json_key)
            
            if not feature_list or not isinstance(feature_list, list):
                raise CmdException('No "%s" found in DSSR output' % json_key)            
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
            elif feature in ('bulges', 'iloops', 'internal', 'junctions', 'sssegments', 'ssSegments', 'multiplets'):
                nts_long = entry.get('nts_long', '')
                if not nts_long:
                    raise CmdException('%s entry missing nts_long field' % feature)
                sel_str = build_selection_from_nts_list([nt.strip() for nt in nts_long.split(',')])
            elif feature == 'nts':
                nt_id = entry.get('nt_id')
                if not nt_id:
                    raise CmdException('Nucleotide entry missing nt_id field')
                chain, resi = parse_nt_id(nt_id)
                sel_str = '(chain %s and resi %s)' % (chain, resi)
            else:
                raise CmdException('Feature type "%s" not yet implemented' % feature)
            cmd.select(name, sel_str, state=state)
            cmd.color('pink', name)
            selected_features.append(feature)
            if not quiet:
                desc = entry.get('name', entry.get('index', index))
                print('dssr_select_v3: created selection "%s" for %s %s (index %d) in state %d' 
                      % (name, feature, desc, index, state))
                print("All features selected: [" + ", ".join(selected_features) + "]" if selected_features else "All features selected: NONE")
    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass
cmd.extend('dssr_select', dssr_select)
cmd.auto_arg[0].update({'dssr_select': cmd.auto_arg[0]['zoom']})