# Bener Dulger
# Date Updated: 01/04/26

'''
This is the plugin for DSSR integration with PyMOL. It allows users to select
a feature in which the feature is then colored to contrast with the rest of the 
structure which is highlighted grey. Listed below are some lines the user could run
into PyMOL's CLI using dssr_select. The tool takes inspiration from the Pseudoknot 
Visualizer. 
'''

#Sample items to run to see current functionality:
#FIRST: Run: fetch 1ehz 
#dssr_select_v3 feature=pairs, index=1, name=pair1
#dssr_select_v3 feature=hairpins, index=1, name=hairpin1, quiet=0
#dssr_select_v3 feature=junctions, index=1, name=junction1, quiet=0
#dssr_select_v3 feature=stems, index=1, name=stem1, quiet=0

#FOR PSEUDOKNOTS (fetch 5TPY first):
#dssr_select_v3 feature=pseudoknot, index=0
#dssr_select_v3 feature=pseudoknot, index=1
#dssr_select_v3 feature=pseudoknot, index=0
# dssr_select_v3 feature=pseudoknot, index=2

from pymol import cmd, CmdException


def unquote(s):  #taken from dssr_block
    s = str(s)
    if s.rstrip()[-1:] not in ('"', "'"):
        return s
    return cmd.safe_eval(s)


def run_dssr_json(pdb_path, exe):
    import subprocess, json

    args = [
        exe,
        '--json',
        '-i=' + pdb_path,
    ]

    try:
        out = subprocess.check_output(args)
    except subprocess.CalledProcessError:
        raise CmdException('"%s" failed' % exe)
    except OSError:
        raise CmdException('Cannot execute exe="%s"' % exe)

    try:
        text = out.decode('utf-8')  #byte to string conversion
        data = json.loads(text)
    except Exception as e:
        raise CmdException('Failed to parse DSSR JSON: %s' % e)
    return data


def parse_nt_id(nt_id):  #sample: A.U2647 - page 38 of manual
    try:
        chain, res = nt_id.split('.', 1)
    except ValueError:
        raise CmdException('Unexpected nt_id format: "%s"' % nt_id)
    digits = []
    for ch in res:
        if ch.isdigit() or ch == '-':
            digits.append(ch)  #only taking in all int values
    if not digits:
        raise CmdException('Could not extract residue number from "%s"' % res)
    resi = ''.join(digits)
    return chain, resi  #generates output ("A",2647)


def parse_atom_id(atom_id):
    try:
        _, rest = atom_id.split('@', 1)
    except ValueError:
        raise CmdException('Unexpected atom_id format: "%s"' % atom_id)
    return parse_nt_id(rest)


def parse_dotbracket_pseudoknots(dotbracket):
    '''    
    Standard pairs: () = Layer 0
    Pseudoknot layers: [] = Layer 1, {} = Layer 2, <> = Layer 3
    '''
    stack_map = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>',
    }
    
    close_to_open = {v: k for k, v in stack_map.items()}
    
    stacks = {} #unmatched open positions
    layers = {} #output being created
    
    layer_assignment = {
        '()': 0,
        '[]': 1,
        '{}': 2,
        '<>': 3
    }
    
    next_layer = 4
    
    #Note: Delete this later
    '''
    Dot-bracket notation on 2TPK:
     (((((((..((((.....[..)))).((((.........)))).....(((((..]....))))))))))))....&..[[[[[.(((((((]]]]].......)))))))..
    '''
    for idx, char in enumerate(dotbracket):
        if char == '.': #skipping over the periods
            continue
        
        if char in stack_map:
            bracket_type = char + stack_map[char]
            if bracket_type not in stacks:
                stacks[bracket_type] = []
            stacks[bracket_type].append(idx)
        
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
            
            if layer not in layers:
                layers[layer] = []
            layers[layer].append((open_idx, idx))
        
        elif char.isalpha():
            if char not in stacks:
                stacks[char] = []
            
            if not stacks[char]:
                stacks[char].append(idx)
            else:
                open_idx = stacks[char].pop()
                
                if char not in layer_assignment:
                    layer_assignment[char] = next_layer
                    next_layer += 1
                
                layer = layer_assignment[char]
                
                if layer not in layers:
                    layers[layer] = []
                layers[layer].append((open_idx, idx))
    return layers


def build_selection_from_layer(layer_pairs, nts_list):
    residues = set()
    for open_idx, close_idx in layer_pairs:
        if open_idx < len(nts_list):
            nt1 = nts_list[open_idx]
            nt1_id = nt1.get('nt_id')
            if nt1_id:
                chain1, resi1 = parse_nt_id(nt1_id)
                residues.add((chain1, resi1))
        
        if close_idx < len(nts_list):
            nt2 = nts_list[close_idx]
            nt2_id = nt2.get('nt_id')
            if nt2_id:
                chain2, resi2 = parse_nt_id(nt2_id)
                residues.add((chain2, resi2))
    
    if not residues:
        return None
    
    clauses = []
    for chain, resi in residues:
        clause = '(chain %s and resi %s)' % (chain, resi)
        clauses.append(clause)
    
    return ' or '.join(clauses)


def build_selection_from_pair(pair_entry):
    nt1 = pair_entry.get('nt1')
    nt2 = pair_entry.get('nt2')
    if not nt1 or not nt2:
        raise CmdException('Pair entry missing nt1 or nt2')
    chain1, resi1 = parse_nt_id(nt1)
    chain2, resi2 = parse_nt_id(nt2)
    residues = {(chain1, resi1), (chain2, resi2)}
    clauses = []
    for chain, resi in residues:
        clause = '(chain %s and resi %s)' % (chain, resi)
        clauses.append(clause)
    
    return ' or '.join(clauses)


def build_selection_from_hbond(hb_entry):
    atom1_id = hb_entry.get('atom1_id')
    atom2_id = hb_entry.get('atom2_id')
    if not atom1_id or not atom2_id:
        raise CmdException('H-bond entry missing atom1_id or atom2_id')
    chain1, resi1 = parse_atom_id(atom1_id)
    chain2, resi2 = parse_atom_id(atom2_id)
    residues = {(chain1, resi1), (chain2, resi2)}
    clauses = []
    for chain, resi in residues:
        clause = '(chain %s and resi %s)' % (chain, resi)
        clauses.append(clause)
    return ' or '.join(clauses)


def build_selection_from_nts_list(nts_list):
    if not nts_list:
        raise CmdException('Empty nucleotide list')
    
    residues = set()
    for nt_id in nts_list:
        chain, resi = parse_nt_id(nt_id)
        residues.add((chain, resi))
    
    clauses = []
    for chain, resi in residues:
        clause = '(chain %s and resi %s)' % (chain, resi)
        clauses.append(clause)
    return ' or '.join(clauses)


def build_selection_from_stem(stem_entry):
    pairs = stem_entry.get('pairs', [])
    if not pairs:
        raise CmdException('Stem has no pairs')
    all_residues = set()
    for pair in pairs:
        nt1 = pair.get('nt1')
        nt2 = pair.get('nt2')
        if nt1:
            chain, resi = parse_nt_id(nt1)
            all_residues.add((chain, resi))
        if nt2:
            chain, resi = parse_nt_id(nt2)
            all_residues.add((chain, resi))
    
    clauses = []
    for chain, resi in all_residues:
        clause = '(chain %s and resi %s)' % (chain, resi)
        clauses.append(clause)
    
    return ' or '.join(clauses)


def build_selection_from_hairpin(hairpin_entry):
    nts_long = hairpin_entry.get('nts_long')
    if not nts_long:
        raise CmdException('Hairpin missing nts_long field')
    
    nts_list = [nt.strip() for nt in nts_long.split(',')]
    return build_selection_from_nts_list(nts_list)


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
    feature = unquote(feature).lower()  #case for "pairs" to pairs and "Pairs" to pairs
    
    feature_map = {
        'pairs': 'pairs',
        'hbonds': 'hbonds',
        'stems': 'stems',
        'helices': 'helices',
        'hairpins': 'hairpins',
        'bulges': 'bulges',
        'iloops': 'iloops',
        'internal': 'iloops',  
        'junctions': 'junctions',
        'ssSegments': 'ssSegments',
        'sssegments': 'ssSegments',  
        'multiplets': 'multiplets',
        'nts': 'nts',
        'pseudoknot': 'pseudoknot',
    }
    
    if feature not in feature_map:
        valid = ', '.join(sorted(feature_map.keys()))
        raise CmdException('Unknown feature "%s". Valid: %s' % (feature, valid))
    
    json_key = feature_map[feature]
    
    if state == 0:
        state = cmd.get_state()
    elif state < 0:
        state = cmd.get_state()
    
    tmpfilepdb = tempfile.mktemp('.pdb')

    try:
        cmd.save(tmpfilepdb, selection, state)
        cmd.color('gray', selection)  #setting base color to grey
        dssr_data = run_dssr_json(tmpfilepdb, exe)
        
        if not quiet:
            print("DSSR top-level JSON keys:")
            print(sorted(dssr_data.keys()))
        
        if feature == 'pseudoknot':
            layer_colors = [
                'blue',      #Standard Watson-Crick pairs
                'pink',       # First pseudoknot layer
                'green',     #  Second pseudoknot layer
                'yellow',    #  Third pseudoknot layer
                'orange',    # Fourth pseudoknot layer
            ]
            
            if 'dbn' not in dssr_data:
                raise CmdException('No dot-bracket notation found in DSSR output')
            
            dbn_data = dssr_data['dbn']
            
            # Extract the actual dot-bracket string from the nested structure
            if isinstance(dbn_data, dict):
                if 'all_chains' in dbn_data and 'sstr' in dbn_data['all_chains']:
                    dotbracket = dbn_data['all_chains']['sstr']
                elif 'sstr' in dbn_data:
                    dotbracket = dbn_data['sstr']
                else:
                    # Searching for sstr dbn in any chain
                    for key, value in dbn_data.items():
                        if isinstance(value, dict) and 'sstr' in value:
                            dotbracket = value['sstr']
                            break
                    else:
                        raise CmdException('Could not find sstr field in dbn data')
            else:
                dotbracket = dbn_data
            
            if not quiet:
                print('Dot-bracket notation:', dotbracket)
            
            if 'nts' not in dssr_data:
                raise CmdException('No nucleotide information found in DSSR output')
            
            nts_list = dssr_data['nts']
            
            layers = parse_dotbracket_pseudoknots(dotbracket)
            
            if not quiet:
                print('Layers found:', sorted(layers.keys()))
                for layer_id in sorted(layers.keys()):
                    print('  Layer %d: %d pairs' % (layer_id, len(layers[layer_id])))
            
            if not layers:
                raise CmdException('No pseudoknot layers found in structure')
            
            if index not in layers:
                available = ', '.join(str(i) for i in sorted(layers.keys()))
                raise CmdException('Layer index %d not found. Available layers: %s' 
                                 % (index, available))
            
            pairs = layers[index]
            color = layer_colors[index % len(layer_colors)]
            
            if show_info:
                print('\n=== Pseudoknot Layer %d ===' % index)
                print('Dot-bracket notation: %s' % dotbracket)
                if index == 0:
                    layer_type = 'Standard Watson-Crick pairs'
                else:
                    layer_type = 'Pseudoknot layer %d' % index
                print('\n%s (%s): %d base pairs' % (layer_type, color, len(pairs)))
            
            sel_str = build_selection_from_layer(pairs, nts_list)
            
            if sel_str is None:
                raise CmdException('Could not build selection for layer %d' % index)
            
            cmd.select(name, sel_str, state=state)
            cmd.color(color, name)
            
            if show_info:
                print('\nBase pairs:')
                for open_idx, close_idx in pairs:
                    if open_idx < len(nts_list) and close_idx < len(nts_list):
                        nt1 = nts_list[open_idx]
                        nt2 = nts_list[close_idx]
                        nt1_id = nt1.get('nt_id', 'unknown')
                        nt2_id = nt2.get('nt_id', 'unknown')
                        print('  %s - %s' % (nt1_id, nt2_id))
                print('\n' + '='*50 + '\n')
            
            if not quiet:
                if index == 0:
                    layer_desc = 'standard pairs'
                else:
                    layer_desc = 'pseudoknot layer %d' % index
                print('dssr_select_v3: created selection "%s" for %s (index %d) with %d base pairs (%s) in state %d' 
                      % (name, layer_desc, index, len(pairs), color, state))
        
        else:
            feature_list = dssr_data.get(json_key, None)
        
            if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
                raise CmdException('No "%s" found in DSSR output' % json_key)
            
            if index < 1 or index > len(feature_list):
                raise CmdException('%s index %d out of range (1..%d)' 
                                 % (feature, index, len(feature_list)))
            
            #conversion system so input of 1 returns first feature: 1 - 1 equals index 0 in python
            #easier for human input
            entry = feature_list[index - 1]
            
            if feature in ('pairs',):
                sel_str = build_selection_from_pair(entry)
            elif feature in ('hbonds',):
                sel_str = build_selection_from_hbond(entry)
            elif feature in ('stems', 'helices'):
                sel_str = build_selection_from_stem(entry)
            elif feature in ('hairpins',):
                sel_str = build_selection_from_hairpin(entry)
            elif feature in ('bulges', 'iloops', 'internal', 'junctions', 'sssegments', 'ssSegments'):
                nts_long = entry.get('nts_long', '')
                if nts_long:
                    nts_list_parsed = [nt.strip() for nt in nts_long.split(',')]
                    sel_str = build_selection_from_nts_list(nts_list_parsed)
                else:
                    raise CmdException('%s entry missing nts_long field' % feature)
            elif feature == 'multiplets':
                nts_long = entry.get('nts_long', '')
                if nts_long:
                    nts_list_parsed = [nt.strip() for nt in nts_long.split(',')]
                    sel_str = build_selection_from_nts_list(nts_list_parsed)
                else:
                    raise CmdException('Multiplet entry missing nts_long field')
            elif feature == 'nts':
                nt_id = entry.get('nt_id')
                if nt_id:
                    chain, resi = parse_nt_id(nt_id)
                    sel_str = '(chain %s and resi %s)' % (chain, resi)
                else:
                    raise CmdException('Nucleotide entry missing nt_id field')
            else:
                raise CmdException('Feature type "%s" not yet implemented' % feature)
            
            cmd.select(name, sel_str, state=state)
            cmd.color('pink', name)  # colors selected feature pink
            selected_features.append(feature)  #adds selected feature to empty list

            if not quiet:
                desc = entry.get('name', entry.get('index', index))
                print('dssr_select_v3: created selection "%s" for %s %s (index %d) in state %d' 
                      % (name, feature, desc, index, state))
                if selected_features:
                    print("All features selected: [" + ", ".join(selected_features) + "]")
                else:
                    print("All features selected: NONE")
    
    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass

cmd.extend('dssr_select_v3', dssr_select_v3)
cmd.auto_arg[0].update({
    'dssr_select_v3': cmd.auto_arg[0]['zoom'],
})