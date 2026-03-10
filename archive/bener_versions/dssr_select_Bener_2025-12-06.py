#Bener Dulger
#Date Updated: 12/6/2025

'''
This is the plugin for DSSR integration with PyMOL. It allows users to select
a feature in which the feature is then colored to contrast with the rest of the 
structure which is highlighted grey. Listed below are some lines the user could run
into PyMOL's CLI using dssr_select. The tool takes inspiration from the Pseudoknot 
Visualizer. 
'''


#Sample items to run to see current functionality:
#FIRST: Run: fetch 1ehz 
#dssr_select feature=pairs, index=1, name=pair1
#dssr_select feature=hairpins, index=1, name=hairpin1, quiet=0
#dssr_select feature=junctions, index=1, name=junction1, quiet=0
#dssr_select feature=hbonds, index=1, name=hbond1, quiet=0
#dssr_select feature=stems, index=1, name=stem1, quiet=0



from pymol import cmd, CmdException


def unquote(s): #taken from dssr_block
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
        text = out.decode('utf-8') #byte to string conversion
        data = json.loads(text)
    except Exception as e:
        raise CmdException('Failed to parse DSSR JSON: %s' % e)
    return data

def parse_nt_id(nt_id): #sample: A.U2647 - page 38 of manual
    try:
        chain, res = nt_id.split('.', 1)
    except ValueError:
        raise CmdException('Unexpected nt_id format: "%s"' % nt_id)
    digits = []
    for ch in res:
        if ch.isdigit() or ch == '-':
            digits.append(ch) #only taking in all int values
    if not digits:
        raise CmdException('Could not extract residue number from "%s"' % res)
    resi = ''.join(digits)
    return chain, resi #generates output ("A",2647)


def parse_atom_id(atom_id):
    try:
        _, rest = atom_id.split('@', 1)
    except ValueError:
        raise CmdException('Unexpected atom_id format: "%s"' % atom_id)
    return parse_nt_id(rest)


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
def dssr_select(selection='all',
                state=-1,
                feature='pairs',
                index=1,
                name='dssr_select',
                exe='x3dna-dssr',
                quiet=1):

    import tempfile, os
    
    state = int(state)
    index = int(index)
    quiet = int(quiet)
    feature = unquote(feature).lower() #case for "pairs" to pairs and "Pairs" to pairs
    
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
        cmd.color('gray', selection) #setting base color to grey
        dssr_data = run_dssr_json(tmpfilepdb, exe)
        feature_list = dssr_data.get(json_key, None)
    
        if feature_list is None or not isinstance(feature_list, list) or len(feature_list) == 0:
            raise CmdException('No "%s" found in DSSR output' % json_key)
        
        if index < 1 or index > len(feature_list):
            raise CmdException('%s index %d out of range (1..%d)' 
                             % (feature, index, len(feature_list)))
        #conversion system so input of 1 returns first feature: 1 - 1 equals index 0 in python
        #easier for human input
        entry = feature_list[index - 1]
        #checking for feature selected by user, calling functions to organize info
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
                nts_list = [nt.strip() for nt in nts_long.split(',')]
                sel_str = build_selection_from_nts_list(nts_list)
            else:
                raise CmdException('%s entry missing nts_long field' % feature)
        elif feature == 'multiplets':
            nts_long = entry.get('nts_long', '')
            if nts_long:
                nts_list = [nt.strip() for nt in nts_long.split(',')]
                sel_str = build_selection_from_nts_list(nts_list)
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
        cmd.color('pink', name) # colors selected feature pink
        selected_features.append(feature) #adds selected feature to empy list

        if not quiet:
            desc = entry.get('name', entry.get('index', index))
            print('dssr_select: created selection "%s" for %s %s (index %d) in state %d' 
                  % (name, feature, desc, index, state))
            if selected_features:
                print("All features selected: [" + ", ".join(selected_features)+ "]")
            else:
                print("All features selected: NONE")
    finally:
        try:
            os.remove(tmpfilepdb)
        except OSError:
            pass

cmd.extend('dssr_select', dssr_select)
cmd.auto_arg[0].update({
    'dssr_select': cmd.auto_arg[0]['zoom'],
})

