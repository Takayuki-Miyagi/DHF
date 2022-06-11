#!/usr/bin/env python3
import sys, os, itertools, subprocess

atoms = ["Ne"]
zetai_list = [1, 2, 3, 4, 5, 6, 7, 8]
#orbitals_list = [f"nonrel-s{n}-p{n}" for n in range(2,3)]
orbitals_list = [f"rel-s{n}-p{n}" for n in range(2,3)]
ARGS = {}
for atom, orbitals, zeta_inv in itertools.product(atoms, orbitals_list, zetai_list):
    ARGS["atom"] = atom
    ARGS["orbitals"] = orbitals
    ARGS["zeta_inv"] = zeta_inv
    if(orbitals.find("nonrel") != -1): ARGS["radial_function_type"] = "NonRel_Laguerre"

    cmd = ' '.join(["./HFTest.exe",] + [f"{x}={ARGS[x]}" for x in ARGS.keys()])
    subprocess.call(cmd, shell=True)    
