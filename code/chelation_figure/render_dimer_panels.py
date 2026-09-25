#!/usr/bin/env python3
"""Render C2 homodimer Zn-cage panels in the manuscript (Fig 5) style.

Reproduces the 4-panel de novo dimer figure. Reads the constructed dimer PDBs
from ./dimers/ (built by build_dimer.py: chain A + chain B = C2 about z, Zn at
the ORIGIN — the PDBs carry no Zn atom, it is added here as a pseudoatom at 0,0,0).

Style: His side-chain C = accent red (#C0392B), imidazole N = blue, Zn = grey
sphere, grey dashed Zn-N coordination bonds, per-protomer cartoon (A blue / B
orange, 30% transparent), white background, ray-traced 300 dpi.

The VIEW for each panel is an explicit sequence of turns from PyMOL's reset frame
(reset = looking down +z, the C2 axis). The locked-in angles below are the ones
chosen for the figure; edit `panels` to re-aim.

RUN (needs open-source PyMOL; e.g. the conda env 'pymolenv'):
    ~/anaconda3/envs/pymolenv/bin/pymol -cq render_dimer_panels.py
Outputs -> ./renders/dimer_<tag>_downz.png
"""
import os
from pymol import cmd

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:      # PyMOL -cq may not set __file__; run from code/chelation_figure/
    HERE = os.getcwd()
DIMERS = os.path.join(HERE, "structures")
OUT    = os.path.join(HERE, "renders"); os.makedirs(OUT, exist_ok=True)

HIS_RED = [0.753, 0.224, 0.169]     # #C0392B  His side-chain carbons
N_BLUE  = [0.145, 0.353, 0.835]     # imidazole nitrogens
ZN_GREY = [0.55, 0.60, 0.70]
CHAIN_A = "0x5b8bd0"                 # protomer A (blue)
CHAIN_B = "0xe08a3c"                 # protomer B (orange)

# tag: (dimer pdb, coordinating (resi, atom) pairs, His resis to show as sticks,
#       [ (axis, degrees), ... ] turn sequence applied after reset+zoom)
PANELS = {
 "zmeg":    ("UA_STRG.15519.1_1_193_357_H32H36_dimer.pdb",  [(32,"ND1"),(36,"NE2")], "32+34+36",
             [("x",-55),("z",10)]),
 "es14265": ("ES_STRG.14265.1_4_293_517_closest_dimer.pdb", [(40,"ND1"),(44,"NE2")], "40+42+44",
             [("x",-55),("z",10)]),
 "dk14275": ("DK_STRG.14275.1_5_1279_1518_closest_dimer.pdb",[(54,"ND1"),(58,"NE2")], "54+56+58",
             [("x",-55),("z",10),("x",60),("y",30)]),      # locked: +60 x then +30 y
 "fi14643": ("FI_STRG.14643.1_8_2494_2775_closest_dimer.pdb",[(5,"ND1"),(7,"NE2")],  "5+7",
             [("x",-55),("z",10),("x",10)]),               # locked: +10 x
}

def render(tag, pdb, coord, his_resis, turns, w=2100, h=1850, dpi=300):
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 1); cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.25); cmd.set("dash_gap", 0.30); cmd.set("dash_width", 3.2)
    cmd.set_color("his_red", HIS_RED); cmd.set_color("n_blue", N_BLUE); cmd.set_color("zn_grey", ZN_GREY)

    cmd.load(os.path.join(DIMERS, pdb), "zd")
    cmd.pseudoatom("zn_obj", pos=[0., 0., 0.], name="ZN", resn="ZN", elem="Zn", chain="Z")

    cmd.hide("everything")
    cmd.show("cartoon", "zd"); cmd.set("cartoon_transparency", 0.30)
    cmd.color(CHAIN_A, "zd and chain A"); cmd.color(CHAIN_B, "zd and chain B")

    hsel = f"zd and resn HIS and resi {his_resis}"
    cmd.show("sticks", f"({hsel}) and not name C+N+O")
    cmd.color("his_red", f"({hsel}) and elem C"); cmd.color("n_blue", f"({hsel}) and elem N")

    cmd.show("spheres", "zn_obj"); cmd.set("sphere_scale", 0.42, "zn_obj"); cmd.color("zn_grey", "zn_obj")
    for ch in ("A", "B"):
        for (resi, atom) in coord:
            nm = f"c_{ch}{resi}"
            cmd.distance(nm, "zn_obj", f"zd and chain {ch} and resi {resi} and name {atom}")
            cmd.hide("labels", nm); cmd.set("dash_color", "grey30", nm)

    cmd.reset(); cmd.zoom("zd or zn_obj", 2)
    for axis, deg in turns:
        cmd.turn(axis, deg)

    outp = os.path.join(OUT, f"dimer_{tag}_downz.png")
    cmd.ray(w, h); cmd.png(outp, dpi=dpi)
    print("wrote", outp)

if __name__ == "__main__":
    for tag, (pdb, coord, his, turns) in PANELS.items():
        render(tag, pdb, coord, his, turns)
