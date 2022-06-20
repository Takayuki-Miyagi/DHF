#!/usr/bin/env python3
import sys, os, itertools, subprocess

atoms = ["Mg"]
#zetai_list = range(2,10)
zetai_list = range(5,6)
orbitals_list = [f"nonrel-s{n}-p{n}" for n in range(3,10)]
#orbitals_list = [f"rel-s{n}-p{n}" for n in range(3,10)]
ARGS = {}
for atom, orbitals, zeta_inv in itertools.product(atoms, orbitals_list, zetai_list):
    ARGS["atom"] = atom
    ARGS["orbitals"] = orbitals
    ARGS["zeta_inv"] = zeta_inv
    ARGS["NMesh"] = 1000
    if(orbitals.find("nonrel") != -1): ARGS["radial_function_type"] = "NonRel_Laguerre"

    cmd = ' '.join(["time","./HFTest.exe",] + [f"{x}={ARGS[x]}" for x in ARGS.keys()])
    subprocess.call(cmd, shell=True)    
