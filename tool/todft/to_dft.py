import numpy as np
from ase import Atoms
import sys
import math

args = sys.argv
if len(args) != 2:
    print("Usage: python to_dft.py <toml file>")
    exit(0)
file_name = args[1]

try:
    from tomli import load as toml_load
except ImportError:
    try: 
        from toml import load as toml_load
        print("WARNING: tomli is not found and toml is found.")
        print("         use of toml package is left for compatibility.")
        print("         please install tomli package.")
        print("HINT: python3 -m pip install tomli")
    except ImportError:
        print("ERROR: tomli is not found")
        print("HINT: python3 -m pip install tomli")
        exit(1)


with open(file_name) as f:
    dirs = toml_load(f)

input_main = dirs["Main"]
input_ase = dirs["ASE"]
input_sol = dirs["Solver"]
pseudopotentials = input_sol.pop("pseudo")

import pprint

print("Parameters for main part")
pprint.pprint(input_main)
print("Parameters for ase part")
pprint.pprint(input_ase)
print("Parameters for solver part")
pprint.pprint(input_sol)
print("Pseudopotentials")
pprint.pprint(pseudopotentials)

try:
    input_h_term = dirs["H_term"]
    h_term_type = 1
    print("Parameters for H_term part")
    pprint.pprint(input_h_term)
except:
    print("Info: No H_term tag")
    h_term_type = 0


### Read/Write file names ###
file_name = input_main["input_xyz_file"]
wfile_name = input_main["output_file_head"]

from ase.io.extxyz import read_xyz, write_xyz

atoms = read_xyz(open(file_name), index=0)
atoms_info = Atoms()
for atom in atoms:
    atoms_info += atom
ratoms = atoms_info.get_positions()
satoms = atoms_info.symbols
natoms = ratoms.shape[0]

z_max = np.max(ratoms, axis=0)[2]
z_min = np.min(ratoms, axis=0)[2]
slab_size = z_max - z_min
print('slab_size=', slab_size)
z_margin = input_main["param"]["z_margin"]
z_atoms = ratoms[:, 2]

z_bottom_most = [z_min - z_margin, z_min + z_margin]
second_min_index = np.where(np.sort(z_atoms) > z_min)[0][0]
z_bottom_2nd_val = np.sort(z_atoms)[second_min_index]
z_bottom_2nd = [z_bottom_2nd_val - z_margin, z_bottom_2nd_val + z_margin]

atoms_bottom_most = Atoms()
for idx in reversed(np.where((z_atoms <= z_bottom_most[1]) & (z_atoms >= z_bottom_most[0]))[0]):
    atoms_bottom_most += atoms_info.pop(idx)
ratoms_bottom_most = atoms_bottom_most.get_positions()
ratoms_bottom_2nd = ratoms[np.where((z_atoms <= z_bottom_2nd[1]) & (z_atoms >= z_bottom_2nd[0]))]

### Add Hydrogen terminators ###
if h_term_type == 1:
    r_SiH = input_h_term["r_SiH"]
    theta = input_h_term["theta"]
    
    dx = r_SiH * math.sin(math.radians(theta * 0.5))
    dy = dx
    dz = r_SiH * math.cos(math.radians(theta * 0.5))
    
    tmp_x_idx = np.arange(ratoms_bottom_most.shape[0])
    tmp_y_idx = np.arange(ratoms_bottom_2nd.shape[0])
    xx, yy = np.meshgrid(tmp_x_idx, tmp_y_idx)
    rxy_atoms_most = ratoms_bottom_most[:, :2]
    rxy_atoms_2nd = ratoms_bottom_2nd[:, :2]
    distances = np.linalg.norm(rxy_atoms_most[xx] - rxy_atoms_2nd[yy], axis=2)
    
    for idx, r2 in enumerate(ratoms_bottom_2nd):
        r1 = ratoms_bottom_most[np.argmin(distances[:, idx])]
        x12 = r1[0] - r2[0]
        y12 = r1[1] - r2[1]
        if (np.abs(x12) < 1e-12 and np.abs(y12) > 1e-12) or (np.abs(x12) > 1e-12 and np.abs(y12) < 1e-12):
            print('(001)surface')
            if np.abs(x12) < 1e-12:
                atoms_info += Atoms(numbers=[1, 1],
                                    positions=[(r2[0], r2[1] - dy, r2[2] - dz), (r2[0], r2[1] + dy, r2[2] - dz)])
            else:
                atoms_info += Atoms(numbers=[1, 1],
                                    positions=[(r2[0] - dx, r2[1], r2[2] - dz), (r2[0] + dx, r2[1], r2[2] - dz)])
        elif np.abs(x12) < 1e-12 and np.abs(y12) < 1e-12:
            print('(111)surface')
            atoms_info += Atoms(numbers=[1], positions=[(r2[0], r2[1], r2[2] - r_SiH)])
else:
    pass   
    
print("Symbols = {}".format(atoms_info.symbols))
print("natoms = {}".format(len(atoms_info)))
### Get cell info ###
#  The input XYZ file is assumed to have the cell info
#  at the last two lines
a1 = input_main["lattice"]["unit_vec"][0]
a2 = input_main["lattice"]["unit_vec"][1]
a3 = np.array([0.0, 0.0, slab_size + input_main["param"]["slab_margin"]])
atoms_info.set_cell(np.array([a1, a2, a3]))
comment = "Lattice= \"{} {} {} {} {} {} {} {} {}\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"".format(a1[0], a1[1],
                                                                                                        a1[2], a2[0],
                                                                                                        a2[1], a2[2],
                                                                                                        a3[0], a3[1],
                                                                                                        a3[2])
write_xyz(open(wfile_name + ".xyz", "w"), atoms_info, comment)
from ase.io import read, write

ase_union = read(wfile_name + ".xyz")
### CIF file ###
ase_union.write(wfile_name + '.cif', format='cif')
### Constraints ###
from ase.constraints import FixAtoms

z_atoms = atoms_info.get_positions()[:, 2]
ase_union.set_constraint(FixAtoms(indices=np.where(z_atoms <= z_bottom_2nd[1])[0]))

### Quantum Espresso ###
from ase.calculators.espresso import Espresso

calculation = input_sol["control"]["calculation"]
kpts = input_ase["kpts"]
command = input_ase["command"]
if calculation == 'scf':
    qe = Espresso(input_data=input_sol,
                  pseudopotentials=pseudopotentials,
                  kpts=kpts,
                  command=command)
elif calculation == 'relax':
    qe = Espresso(input_data=input_sol,
                  pseudopotentials=pseudopotentials,
                  kpts=kpts,
                  command=command)
elif calculation == 'bands':
    qe = Espresso(input_data=input_sol,
                  pseudopotentials=pseudopotentials,
                  kpts=kpts,
                  command=command)
ase_union.calc = qe
ase_union.calc.write_input(atoms_info)
try:
    ase_union.get_potential_energy()  # run pw.x
except:
    print("Calculation of get_potential_energy is not normally finished.")
    pass
