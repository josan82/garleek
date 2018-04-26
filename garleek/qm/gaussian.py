#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Garleek - Gaussian bridge
"""

from __future__ import print_function, absolute_import, division
from collections import OrderedDict
import numpy as np


def patch_gaussian_input(filename, atom_types, engine='tinker', forcefield=None):

    def _is_route(line):
        return line.startswith('#') and 'external=' in line.lower()

    def _patch_mm_keyword(line):
        command = 'garleek-backend --qm gaussian --mm {}'.format(engine)
        if forcefield:
            command += ' --ff {}'.format(forcefield)
        return line.replace('garleek', '"{}"'.format(command))

    def _patch_atom_type(line):
        fields = line.split()
        atom_fields = fields[0].split('-')
        atom_type = atom_fields[1]
        line = line.replace(atom_type, atom_types[atom_type], 1)
        if len(fields) > 6:
            link_atom = fields[6]
            link_atom_type = link_atom.split('-')[1]
            line = line.replace(link_atom_type, atom_types[link_atom_type])
        return line

    skipped_mult_charges = False
    section = 0
    lines = []
    with open(filename) as f:
        for line in f:
            orig_line, line = line, line.strip()
            if line.startswith('!'):
                pass
            elif not line:
                section += 1
            elif _is_route(line):
                orig_line = _patch_mm_keyword(orig_line)
            elif line and section == 2:
                if skipped_mult_charges:
                    orig_line = _patch_atom_type(orig_line)
                else:
                    skipped_mult_charges = True
            lines.append(orig_line)
    return ''.join(lines)


def parse_gaussian_EIn(ein_filename):
    """
    Parse the `*.EIn`file produced by Gaussian `external` keyword.

    This file contains the following data (taken from http://gaussian.com/external)

    n_atoms  derivatives-requested  charge  spin
    atom_name  x  y  z  MM-charge [atom_type]
    atom_name  x  y  z  MM-charge [atom_type]
    atom_name  x  y  z  MM-charge [atom_type]

    ...

    `derivatives-requested` can be 0 (energy only), 1 (first derivatives)
    or 2 (second derivatives).
    """
    with open(ein_filename) as f:
        n_atoms, derivatives, charge, spin = list(map(int, next(f).split()))
        atoms = OrderedDict()
        for i in range(n_atoms):
            fields = next(f).strip().split()
            atom_element = fields[0]
            atom_type = fields[5] if len(fields) == 6 else None
            x, y, z, mm_charge = list(map(float, fields[1:5]))
            atoms[i+1] = {'element': 'E'+atom_element,
                          'type': atom_type,
                          'xyz': np.array([x, y, z]),
                          'mm_charge': mm_charge}

        line = next(f)  # Skip the "connectivity" header
        if 'connectivity' in line.strip().lower():
            line = next(f)
        bonds = OrderedDict()
        while line.strip():
            fields = line.strip().split()
            bonds[int(fields[0])] = bond_list = []
            for to_atom, bond_index in zip(fields[1::2], fields[2::2]):
                bond_list.append((int(to_atom), float(bond_index)))
            line = next(f, '')

    return {'n_atoms': n_atoms,
            'derivatives': derivatives,
            'charge': charge,
            'spin': spin,
            'atoms': atoms,
            'bonds': bonds}


def prepare_gaussian_EOu(n_atoms, energy, dipole_moment, gradients=None, hessian=None,
                         polarizability=None, dipole_polarizability=None):
    """
    Generate the `*.EOu` file Gaussian expects after `external` launch.

    After performing the MM calculations, Gaussian expects a file with the
    following information (all in atomic units; taken from
    http://gaussian.com/external)

    Items                        Pseudo Code                            Line Format
    -------------------------------------------------------------------------------
    energy, dipole-moment (xyz)  E, Dip(I), I=1,3                       4D20.12
    gradient on atom (xyz)       FX(J,I), J=1,3; I=1,NAtoms             3D20.12
    polarizability               Polar(I), I=1,6                        3D20.12
    dipole derivatives           DDip(I), I=1,9*NAtoms                  3D20.12
    force constants              FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2  3D20.12
    """
    lines = [[energy] + list(dipole_moment)]
    template = '{: 20.12e}'
    if gradients is not None:
        for gradient in gradients:
            lines.append(gradient)
    if hessian is not None:
        if polarizability is None:
            polarizability = np.zeros(6)
        if dipole_polarizability is None:
            dipole_polarizability = np.zeros(9*n_atoms)
        for i in range(0, polarizability.size, 3):
            lines.append(polarizability[i:i+3])
        for i in range(0, dipole_polarizability.size, 3):
            lines.append(dipole_polarizability[i:i+3])

        for i in range(0, hessian.size, 3):
            lines.append(hessian[i:i+3])
    lines.append([])  # Gaussian is very peculiar about blank lines
    return '\n'.join([(template*len(line)).format(*line) for line in lines])
