#!/usr/bin/env python

from __future__ import print_function, division, absolute_import
import os
import sys
import shutil
import re
import pytest
import cclib
from contextlib import contextmanager
from distutils.spawn import find_executable
from tempfile import mkdtemp
from subprocess import call
import numpy as np

from garleek.cli import frontend_app as frontend_garleek, _extant_file_types, _extant_file_prm
here = os.path.abspath(os.path.dirname(__file__))
data = os.path.join(here, 'data')
gaussian_exe = find_executable('g16') or find_executable('g09') or 'g16'
WORKING_DIR = os.getcwd()


def isclose(a, b, rel_tol=1e-06, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def oniom_energy(path):
    energy = None
    with open(path) as f:
        for line in f:
            if 'extrapolated' in line:
                energy = float(line.split()[-1])
    return energy


def check_errors(path):
    last_lines = []
    with open(path) as f:
        for line in f:
            last_lines.append(line)
            last_lines = last_lines[-6:]
            if 'Failed to open output file from external program' in line:
                for last_line in last_lines[:-1]:
                    title, g_output, _ = last_line.split('"')
                    if not os.path.isfile(g_output):
                        continue
                    g_output_copy = os.path.basename(path) + os.path.splitext(g_output)[1]
                    shutil.copy(g_output,  os.path.join(WORKING_DIR, 'outputs', g_output_copy))
                return False
    return True


@contextmanager
def temporary_directory(enter=True, remove=True, **kwargs):
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = mkdtemp(prefix='garleek', **kwargs)
    if enter:
        cwd = os.getcwd()
        os.chdir(temp_dir)
    yield temp_dir
    if enter:
        os.chdir(cwd)
    if remove:
        shutil.rmtree(temp_dir)


@pytest.mark.parametrize("directory", sorted(next(os.walk(data))[1]))
def test_gaussian_tinker(directory):
    if directory.endswith('UFF'):
        pytest.skip()

    with temporary_directory(remove=False) as tmp:
        # Original data
        data_original = os.path.join(data, directory)
        outfile_original = os.path.join(data_original, directory + '.out')
        infile_original = os.path.join(data_original, directory + '.in')
        # Copied tmp paths
        data_copy = os.path.join(tmp, directory)
        infile_copy = os.path.join(data_copy, directory + '.in')
        shutil.copytree(data_original, data_copy)
        os.chdir(data_copy)
        # We are now in /tmp/garleek*****/A_1 or similar, which
        # contains copies of the original Gaussian files and types
        # guess forcefield specified in atom.types
        ff = 'mm3.prm'
        with open('atom.types') as f:
            for line in f:
                if line.startswith('# forcefield:'):
                    ff = line.split(':', 1)[1].strip()
        # patch inputfile
        garleek_in = frontend_garleek(infile_copy, qm='gaussian', mm='tinker',
                                      ff=_extant_file_prm(ff), types='atom.types')

        call([gaussian_exe, garleek_in])
        garleek_out = os.path.splitext(garleek_in)[0] + '.log'

        # Save output in working dir
        assert os.path.isfile(garleek_out)
        if not os.path.exists(os.path.join(WORKING_DIR, 'outputs')):
            os.mkdir(os.path.join(WORKING_DIR, 'outputs'))
        shutil.copy(garleek_out, os.path.join(WORKING_DIR, 'outputs', os.path.basename(garleek_out)))
        assert check_errors(garleek_out)
        # Check values
        cc_original = cclib.ccopen(outfile_original).parse()
        cc_calculated = cclib.ccopen(garleek_out).parse()
        assert isclose(cc_original.scfenergies[-1], cc_calculated.scfenergies[-1])
        assert np.sqrt(np.mean(np.square(cc_original.atomcoords[-1]-cc_calculated.atomcoords[-1]))) < 0.001
        assert isclose(oniom_energy(outfile_original), oniom_energy(garleek_out))