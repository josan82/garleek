import re
import os
import sys

atom_rx = r'^atom\s+([0-9]+)\s+([0-9]+)\s+(\w+)\s+"(.*)"\s+([0-9]+)\s+([0-9.]+)\s+([0-9]+)$'

AA = {
    'alanine': 'ALA',
    'arginine': 'ARG',
    'asparagine': 'ASN',
    'aspartic': 'ASP',
    'cysteine': 'CYS',
    'glutamic': 'GLU',
    'glutamine': 'GLN',
    'glycine': 'GLY',
    'histidine': 'HIS',
    'isoleucine': 'ILE',
    'leucine': 'LEU',
    'lysine': 'LYS',
    'methionine': 'MET',
    'phenylalanine': 'PHE',
    'proline': 'PRO',
    'serine': 'SER',
    'threonine': 'THR',
    'tryptophan': 'TRP',
    'tyrosine': 'TYR',
    'valine': 'VAL',
}

d = {}
with open(sys.argv[1]) as f:
    for line in f:
        search = re.search(atom_rx, line)
        if search:
            atype, aclass, name, descr, symbol, mass, valence = search.groups()
            residue = descr.split()[0].lower()
            if residue in AA:
                restype = AA[residue]
                d[restype + '_' + name] = atype

with open('atom.types', 'w') as f:
    f.write('\n'.join(k + ' ' + v for (k,v) in sorted(d.items())))