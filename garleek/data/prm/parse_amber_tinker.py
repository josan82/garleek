import re
import os
import sys

atom_rx = r'^atom\s+([0-9]+)\s+([0-9]+)\s+(\w+)\s+"(.*)"\s+([0-9]+)\s+([0-9.]+)\s+([0-9]+)$'

AA = {
    'alanine':              'ALA',
    'arginine':             'ARG',
    'asparagine':           'ASN',
    'aspartic acid':        'ASP',
    'cysteine':             'CYS',
    'cysteine (-sh)':       'CYS', # Review
    'cystine (-ss-)':       'CYX',
    'glutamic acid':        'GLU',
    'glutamine':            'GLN',
    'glycine':              'GLY',
    'histidine':            'HIS',
    'histidine (+)':        'HIP',
    'histidine (hd)':       'HID',
    'histidine (he)':       'HIE',
    'isoleucine':           'ILE',
    'leucine':              'LEU',
    'lysine':               'LYS',
    'methionine':           'MET',
    'phenylalanine':        'PHE',
    'proline':              'PRO',
    'serine':               'SER',
    'threonine':            'THR',
    'tryptophan':           'TRP',
    'tyrosine':             'TYR',
    'valine':               'VAL',

    'ornithine':            'ORN',
    'methylalanine':        '',
    'pyroglutamate':        '',
    'formyl':               '',
    "acetyl":               'ACE',
    "c-term amide":         'NHE',
    "n-meamide":            'NME',

    "n-term aib":           'NAIB',
    "n-term ala":           'NALA',
    "n-term arg":           'NARG',
    "n-term asn":           'NASN',
    "n-term asp":           'NASP',
    "n-term cys (-sh)":     'NCYS',
    "n-term cys (-ss-)":    'NCYS',
    "n-term gln":           'NGLN',
    "n-term glu":           'NGLU',
    "n-term gly":           'NGL',
    "n-term his (+)":       'NHIS',
    "n-term his (hd)":      'NHIS',
    "n-term his (he)":      'NHIS',
    "n-term ile":           'NILE',
    "n-term leu":           'NLEU',
    "n-term lys":           'NLYS',
    "n-term met":           'NMET',
    "n-term orn":           'NORN',
    "n-term phe":           'NPH',
    "n-term pro":           'NPRO',
    "n-term ser":           'NSER',
    "n-term thr":           'NTHR',
    "n-term trp":           'NTRP',
    "n-term tyr":           'NTYR',
    "n-term val":           'NVAL',

    "c-term aib":           'CAIB',
    "c-term ala":           'CALA',
    "c-term arg":           'CARG',
    "c-term asn":           'CASN',
    "c-term asp":           'CAS',
    "c-term cys (-sh)":     'CCYS',
    "c-term cys (-ss-)":    'CCYX',
    "c-term gln":           'CGLN',
    "c-term glu":           'CGLU',
    "c-term gly":           'CGL',
    "c-term his (+)":       'CHIS',
    "c-term his (hd)":      'CHIS',
    "c-term his (he)":      'CHIS',
    "c-term ile":           'CILE',
    "c-term leu":           'CLEU',
    "c-term lys":           'CLYS',
    "c-term met":           'CMET',
    "c-term orn":           'CORN',
    "c-term phe":           'CPHE',
    "c-term pro":           'CPRO',
    "c-term ser":           'CSER',
    "c-term thr":           'CTHR',
    "c-term trp":           'CTRP',
    "c-term tyr":           'CTYR',
    "c-term val":           'CVAL',

    "r-adenosine":          '',
    "r-guanosine":          '',
    "r-cytosine":           '',
    "r-uracil":             '',
    "d-adenosine":          '',
    "d-guanosine":          '',
    "d-cytosine":           '',
    "d-thymine":            '',
    "r-phosphodiester":     '',
    "r-5'-hydroxyl":        '',
    "r-5'-phosphate":       '',
    "r-3'-hydroxyl":        '',
    "r-3'-phosphate":       '',
    "d-phosphodiester":     '',
    "d-5'-hydroxyl":        '',
    "d-5'-phosphate":       '',
    "d-3'-hydroxyl":        '',
    "d-3'-phosphate":       '',
    "tip3p":                'HOH',
    "li+":                  'LI',
    "na+":                  'NA',
    "k+":                   'K',
    "rb+":                  'RB',
    "cs+":                  'CS',
    "mg+2":                 'MG',
    "ca+2":                 'CA',
    "zn+2":                 'ZN',
    "ba+2":                 'BA',
    "cl-":                  'CL',
}
ALIASES = {
    'CT': ('CX', 'CA'),
    'H1': ('HP', 'HC', 'HA'),
    'O2': ('O',),
}
# default assignments
d = { 
    'H1': '12', # link atom
    '1':  '12',
    'HC': '14',
}
with open(sys.argv[1]) as f:
    for line in f:
        search = re.search(atom_rx, line)
        if search:
            atype, aclass, name, descr, symbol, mass, valence = search.groups()
            residue = ' '.join(descr.split()[:-1]).lower()
            if residue in AA:
                restype = AA[residue].upper()
                if restype:
                    d[restype + '_' + name] = atype
                    for alias in ALIASES.get(name, []):
                        d[restype + '_' + alias] = atype


try:
    output = sys.argv[2]
except IndexError:
    output = 'atom.types'
with open(output, 'w') as f:
    f.write('\n'.join(k + ' ' + v for (k,v) in sorted(d.items())))