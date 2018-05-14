import re
import os
import sys


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
    'OE': ('OE1', 'OE2'),
    'OD': ('OD1', 'OD2'),
    'NH': ('NH1', 'NH2'),
    'NE': ('NE2',),
    'ND': ('ND1',),
    'HA': ('HA1', 'HA2', 'HA3'),
    'HB': ('HB1', 'HB2', 'HB3'),
    'HC': ('HC1', 'HC2', 'HC3'),
    'HD': ('HD1', 'HD2', 'HD3'),
    'HD1': ('HD11', 'HD12', 'HD13'),
    'HD2': ('HD21', 'HD22', 'HD23'),
    'HE': ('HE1', 'HE2', 'HE3'),
    'HE1': ('HE11', 'HE12'),
    'HE2': ('HE21', 'HE22'),
    'HG': ('HG1', 'HG2', 'HG3',),
    'HG1': ('HG11', 'HG12', 'HG13'),
    'HG2': ('HG21', 'HG22', 'HG23',),
    'HH': ('HH1', 'HH2', 'HH11', 'HH12', 'HH21', 'HH22'),
    'HH1': ('HH11', 'HH12'),
    'HH2': ('HH21', 'HH22'),
    'HZ': ('HZ1', 'HZ2', 'HZ3'),
    'CE': ('CE1', 'CE2'),
    'CD': ('CD1', 'CD2'),
    'HN': ('H'),
    'OXT': ('O'),
}

def main(prm_file):
    biotype_rx = r'^biotype\s+\d+\s+(\S+)\s+"(.*)"\s+(\d+)'
    # default assignments
    d = {
        'H1': '12',  # link atom
        '1':  '12',
        'HC': '14',
    }
    with open(prm_file) as f:
        for line in f:
            search = re.search(biotype_rx, line)
            if search:
                atom_type, residue_str, tinker_type = search.groups()
                if residue_str in AA:
                    restype = AA[residue_str.lower()].upper()
                    if restype:
                        d[restype + '_' + atom_type] = tinker_type
                        for alias in ALIASES.get(atom_type, ''):
                            d[restype + '_' + alias] = tinker_type
    return d

if __name__ == '__main__':
    try:
        d = main(sys.argv[1])
    except IndexError:
        sys.exit('Usage: biotyper.py forcefield.prm output.types')
    else:
        print('\n'.join('{} {}'.format(k, v) for (k,v) in sorted(d.items())))
