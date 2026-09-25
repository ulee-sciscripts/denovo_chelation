#!/usr/bin/env python3
"""Render the FI14643 His2Cys2 Zn-cage panel from the Boltz model (model_3),
matching the manuscript style of render_dimer_panels.py.

Site (verified): inter-chain His2Cys2 = Cys62(SG) + His89(NE2) from EACH chain,
one Zn (Metal3D 0.9993). Zn coords come from the model; this is a low-confidence
predicted capacity (iPTM 0.12) rendered for illustration.

RUN (needs open-source PyMOL):
    ~/anaconda3/envs/pymolenv/bin/pymol -cq render_fi14643_his2cys2.py
Out -> ./renders/dimer_fi14643_his2cys2_downz.png
"""
import os
from pymol import cmd

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:      # PyMOL -cq may not set __file__; run from code/chelation_figure/
    HERE = os.getcwd()
CIF  = os.path.join(HERE, "structures", "FI_STRG_14643_1_8_2494_2775_bestZn_model3.cif")
OUT  = os.path.join(HERE, "renders"); os.makedirs(OUT, exist_ok=True)

HIS_RED = [0.753, 0.224, 0.169]   # #C0392B side-chain carbons
N_BLUE  = [0.145, 0.353, 0.835]   # imidazole N
S_YELL  = [0.95, 0.78, 0.18]      # Cys sulfur
ZN_GREY = [0.55, 0.60, 0.70]
CHAIN_A = "0x5b8bd0"
CHAIN_B = "0xe08a3c"

def main():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 1); cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.25); cmd.set("dash_gap", 0.30); cmd.set("dash_width", 3.2)
    cmd.set_color("his_red", HIS_RED); cmd.set_color("n_blue", N_BLUE)
    cmd.set_color("s_yell", S_YELL); cmd.set_color("zn_grey", ZN_GREY)

    cmd.load(CIF, "m")
    cmd.remove("solvent")
    chains = cmd.get_chains("m and polymer")
    print("chains:", chains)

    cmd.hide("everything")
    cmd.show("cartoon", "m and polymer"); cmd.set("cartoon_transparency", 0.30)
    ca, cb = (chains + ["A", "B"])[:2]
    cmd.color(CHAIN_A, f"m and chain {ca}")
    cmd.color(CHAIN_B, f"m and chain {cb}")

    site = "m and ((resi 62 and resn CYS) or (resi 89 and resn HIS))"
    cmd.show("sticks", f"({site}) and not name C+N+O")
    cmd.color("his_red", f"({site}) and elem C")
    cmd.color("n_blue",  f"({site}) and elem N")
    cmd.color("s_yell",  f"({site}) and elem S")

    # Zn sphere at its modeled position
    cmd.show("spheres", "m and elem Zn")
    cmd.set("sphere_scale", 0.45, "m and elem Zn")
    cmd.color("zn_grey", "m and elem Zn")

    # dashed coordination bonds Zn -> each donor
    znsel = "m and elem Zn"
    for i, dsel in enumerate([
        f"m and chain {ca} and resi 62 and name SG",
        f"m and chain {cb} and resi 62 and name SG",
        f"m and chain {ca} and resi 89 and name NE2",
        f"m and chain {cb} and resi 89 and name NE2",
    ]):
        cmd.distance(f"coord{i}", znsel, dsel)
        cmd.color("gray50", f"coord{i}"); cmd.hide("labels", f"coord{i}")

    # view: center the Zn site, show both protomers
    cmd.orient(f"({site}) or ({znsel})")
    cmd.zoom("m and polymer", 6)      # a hair more padding than before (was 3)
    cmd.turn("y", 15)

    cmd.set("ray_trace_mode", 0)
    cmd.ray(2100, 1850)
    out = os.path.join(OUT, "dimer_fi14643_his2cys2_downz.png")
    cmd.png(out, dpi=300)
    print("wrote", out)

main()
