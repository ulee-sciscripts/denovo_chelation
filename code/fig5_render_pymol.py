#!/usr/bin/env python3
"""Render all ray-traced PyMOL panels for Figure 5.
Run headless:  pymol -cq fig5_render_pymol.py

Unified austere scheme:
  His side-chain carbons = accent RED; imidazole N = blue; metal = grey sphere.
  Panel-A fold cartoon coloured by pLDDT (AF/ESMFold bands); site zooms / cluster
  panels use a faint grey cartoon for context. White background, 300 dpi.

Reads the input structures from $FIG5_SCR (default ./fig5_render):
  cand3.pdb                        (ZMEG apo ESMFold model)
  calprotectin_heterodimer.pdb     (S100A8/S100A9)

Produces (in the same dir), consumed by fig5_metalsite.py:
  zmeg_overview.png / zmeg_overview_markers.png       (Panel A: fold, pLDDT, His run)
  zmeg_site_ion_nolabel.png                           (Panel A: His32-His36 + modeled Zn;
                                                        His labels added off-structure by the compose step)
  zmeg_cluster.png / zmeg_cluster_markers.png         (Panel C: ZMEG His32/His36 cage, apo)
  calpro_cluster.png / calpro_cluster_markers.png     (Panel C: calprotectin site-2 His cage)
"""
import os
from pymol import cmd

SCR = os.environ.get("FIG5_SCR", "fig5_render")   # scratch dir: input PDBs (cand3.pdb, calprotectin_heterodimer.pdb) + rendered panels

AF_COLORS = {
    "af_vlow": (1.000, 0.490, 0.271),  # <50  very low
    "af_low":  (1.000, 0.859, 0.075),  # 50-70 low
    "af_conf": (0.396, 0.796, 0.953),  # 70-90 confident
    "af_high": (0.000, 0.325, 0.839),  # >90  very high
}
HIS_RED = (0.753, 0.224, 0.169)   # #C0392B accent
N_BLUE  = (0.145, 0.353, 0.835)   # imidazole nitrogen
ZN_GREY = (0.55, 0.60, 0.70)

W, H = 2000, 1800

# marker-pass pure colours (index-matched to residue lists), must match compose script
MARKER_RGB = [
    (1.0, 0.0, 0.0), (0.0, 0.7, 0.0), (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0), (0.0, 0.85, 0.85), (1.0, 0.6, 0.0),
]


def global_settings():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.set("hash_max", 300)
    cmd.set("ambient", 0.42)
    cmd.set("direct", 0.55)
    cmd.set("specular", 0.25)
    cmd.set("spec_power", 60)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_highlight_color", "grey60")
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("stick_radius", 0.22)
    cmd.set("label_size", 26)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_font_id", 7)
    cmd.set("dash_gap", 0.35)
    cmd.set("dash_width", 3.0)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)
    for name, rgb in AF_COLORS.items():
        cmd.set_color(name, list(rgb))
    cmd.set_color("his_red", list(HIS_RED))
    cmd.set_color("n_blue", list(N_BLUE))
    cmd.set_color("zn_grey", list(ZN_GREY))


def color_by_plddt(obj):
    cmd.alter(obj, "b = b * 100.0")
    cmd.rebuild()
    cmd.color("af_vlow", f"{obj}")
    cmd.color("af_low",  f"{obj} and b > 50")
    cmd.color("af_conf", f"{obj} and b > 70")
    cmd.color("af_high", f"{obj} and b > 90")


def style_his(sel):
    """His side chains: red carbons, blue N."""
    cmd.show("sticks", f"({sel}) and (sidechain or name CA+CB)")
    cmd.set("stick_color", "his_red", f"({sel}) and elem C")
    cmd.color("n_blue", f"({sel}) and elem N")
    cmd.color("red", f"({sel}) and elem O")


def marker_pass(obj, res_list, marker_png):
    """Flat pure-colour spheres at each residue CB (view-independent label anchors).
    res_list = list of (chain, resi)."""
    view = cmd.get_view()
    cmd.hide("everything", obj)
    for i, (ch, r) in enumerate(res_list):
        anchor = f"{obj} and chain {ch} and resi {r} and name CB"
        if cmd.count_atoms(anchor) == 0:
            anchor = f"{obj} and chain {ch} and resi {r} and name CA"
        cname = f"mk_{i}"
        cmd.set_color(cname, list(MARKER_RGB[i % len(MARKER_RGB)]))
        cmd.show("spheres", anchor)
        cmd.set("sphere_scale", 0.75, anchor)
        cmd.color(cname, anchor)
    cmd.set("ambient", 1.0); cmd.set("direct", 0.0); cmd.set("specular", 0.0)
    cmd.set("ray_shadows", 0); cmd.set("antialias", 0)
    cmd.set_view(view)
    cmd.png(marker_png, width=W, height=H, dpi=300, ray=1)
    print("wrote", marker_png)
    cmd.set("ambient", 0.42); cmd.set("direct", 0.55); cmd.set("specular", 0.25)
    cmd.set("ray_shadows", 1); cmd.set("antialias", 2)


# ---------------- Panel A: ZMEG fold overview ----------------
def zmeg_overview():
    cmd.reinitialize(); global_settings()
    cmd.load(f"{SCR}/cand3.pdb", "ZMEG"); cmd.remove("solvent")
    his = [27, 29, 32, 34, 36]
    cmd.hide("everything", "ZMEG"); cmd.show("cartoon", "ZMEG")
    color_by_plddt("ZMEG")
    hsel = f"ZMEG and resi {'+'.join(map(str, his))}"
    cmd.set("stick_radius", 0.26)
    style_his(hsel)
    cmd.orient("ZMEG"); cmd.zoom("ZMEG", buffer=6.0); cmd.turn("y", 8)
    out = f"{SCR}/zmeg_overview.png"
    cmd.png(out, width=W, height=H, dpi=300, ray=1); print("wrote", out)
    marker_pass("ZMEG", [("A", r) for r in his], f"{SCR}/zmeg_overview_markers.png")


# ---------------- Panel A: ZMEG site + modeled ion ----------------
def zmeg_site_ion():
    cmd.reinitialize(); global_settings()
    cmd.load(f"{SCR}/cand3.pdb", "ZMEG"); cmd.remove("solvent")
    r1, r2 = 32, 36
    ion = (10.078, 7.988, 12.370)
    n1 = "ZMEG and resi 32 and name ND1"
    n2 = "ZMEG and resi 36 and name NE2"
    site = f"ZMEG and resi {r1}+{r2}"
    cmd.hide("everything", "ZMEG")
    cmd.show("cartoon", "ZMEG"); cmd.set("cartoon_transparency", 0.55, "ZMEG")
    cmd.color("grey80", "ZMEG")
    style_his(site)
    cmd.show("spheres", f"({n1}) or ({n2})")
    cmd.set("sphere_scale", 0.30, f"({n1}) or ({n2})")
    cmd.color("n_blue", f"({n1}) or ({n2})")
    cmd.distance("nnd", n1, n2); cmd.hide("labels", "nnd")
    cmd.set("dash_color", "grey30", "nnd")
    cmd.pseudoatom("ion_ps", pos=list(ion))
    cmd.show("spheres", "ion_ps"); cmd.set("sphere_scale", 0.55, "ion_ps")
    cmd.color("zn_grey", "ion_ps")
    for nm, nsel in (("c1", n1), ("c2", n2)):
        cmd.distance(nm, "ion_ps", nsel); cmd.hide("labels", nm)
        cmd.set("dash_color", "orange", nm); cmd.set("dash_gap", 0.30, nm)
        cmd.set("dash_width", 3.5, nm)
    # NB: His32/His36 are NOT baked in here -- they are added off-structure with
    # leader lines by fig5_metalsite.py (the baked labels overlapped the ribbon).
    focus = "((ZMEG and resi 32+36) or ion_ps)"
    cmd.orient(focus); cmd.zoom(focus, buffer=2.2)
    out = f"{SCR}/zmeg_site_ion_nolabel.png"
    cmd.png(out, width=W, height=H, dpi=300, ray=1); print("wrote", out)


# ---------------- Panel C clusters ----------------
def render_cluster(pdb, obj, res_list, out_prefix, zoom_buffer=2.6,
                   mark_N=True):
    """Austere His-cage render for the comparison panel. res_list=[(chain,resi),...]"""
    cmd.reinitialize(); global_settings()
    cmd.load(pdb, obj); cmd.remove("solvent")
    cmd.hide("everything", obj)
    sel_parts = [f"(chain {c} and resi {r})" for c, r in res_list]
    site = f"{obj} and (" + " or ".join(sel_parts) + ")"
    # faint local cartoon context
    cmd.show("cartoon", obj)
    cmd.set("cartoon_transparency", 0.72, obj)
    cmd.color("grey80", obj)
    style_his(site)
    if mark_N:
        nsel = f"({site}) and name ND1+NE2"
        cmd.show("spheres", nsel)
        cmd.set("sphere_scale", 0.26, nsel)
        cmd.color("n_blue", nsel)
    cmd.orient(site); cmd.zoom(site, buffer=zoom_buffer)
    out = f"{SCR}/{out_prefix}.png"
    cmd.png(out, width=W, height=H, dpi=300, ray=1); print("wrote", out)
    marker_pass(obj, res_list, f"{SCR}/{out_prefix}_markers.png")


def zmeg_cluster():
    render_cluster(f"{SCR}/cand3.pdb", "ZMEGC", [("A", 32), ("A", 36)],
                   "zmeg_cluster", zoom_buffer=2.4)


def calpro_cluster():
    render_cluster(f"{SCR}/calprotectin_heterodimer.pdb", "CAL",
                   [("A", 17), ("A", 27), ("B", 91), ("B", 95), ("B", 103)],
                   "calpro_cluster", zoom_buffer=2.2)


if __name__ == "__main__" or True:
    zmeg_overview()
    zmeg_site_ion()
    zmeg_cluster()
    calpro_cluster()
    print("ALL FIG5 PANELS DONE")
