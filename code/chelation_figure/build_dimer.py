#!/usr/bin/env python3
"""Explicit geometric construction of a symmetric (C2) homodimer Zn site
for de novo His-motif peptides. Idealized-symmetry cross-check.

For each target:
  - pick coordinating His pair (specified or closest imidazole-N pair)
  - build ideal tetrahedral Zn template (Zn origin, 4 N at tet vertices x 2.05 A)
  - rigid-body place monomer so its 2 coordinating N sit on 2 template vertices
    (leaves 1 free rotational DOF about the N-N axis)
  - chain B = C2 (about z) of chain A -> fills the other 2 vertices
  - jigger scan: rotate chain A about N-N axis 0..360, min inter-chain clash
  - report realized coordination geometry; write dimer PDB (no Zn) for Metal3D
"""
import sys, os, itertools
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from scipy.spatial import cKDTree

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    HERE = os.getcwd()
# Monomer ESMFold models (external; not committed). Set MONOMER_PDBDIR, or place them under structures/monomers/.
PDBDIR = os.environ.get("MONOMER_PDBDIR", os.path.join(HERE, "structures", "monomers"))
OUT = os.path.join(HERE, "structures")
os.makedirs(OUT, exist_ok=True)

ZN_N = 2.05  # Zn-N distance
# tetrahedral unit vectors
TET = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float)
TET = TET / np.linalg.norm(TET, axis=1, keepdims=True) * ZN_N
VERT_SPACING = np.linalg.norm(TET[0]-TET[1])  # ~3.348
# C2 about z maps v1,v2 -> v4,v3 (occupied -> empty)
C2 = np.diag([-1.0, -1.0, 1.0])

def rot_align(u, v):
    """Rotation matrix aligning unit vector u onto unit vector v."""
    u = u/np.linalg.norm(u); v = v/np.linalg.norm(v)
    c = np.dot(u, v)
    if c > 0.9999999:
        return np.eye(3)
    if c < -0.9999999:
        # 180 deg: pick any axis perpendicular to u
        a = np.array([1,0,0]) if abs(u[0])<0.9 else np.array([0,1,0])
        axis = np.cross(u, a); axis/=np.linalg.norm(axis)
        return rot_axis(axis, np.pi)
    axis = np.cross(u, v); s = np.linalg.norm(axis); axis/=s
    ang = np.arccos(np.clip(c,-1,1))
    return rot_axis(axis, ang)

def rot_axis(axis, ang):
    axis = axis/np.linalg.norm(axis)
    x,y,z = axis; c=np.cos(ang); s=np.sin(ang); C=1-c
    return np.array([
        [c+x*x*C,   x*y*C-z*s, x*z*C+y*s],
        [y*x*C+z*s, c+y*y*C,   y*z*C-x*s],
        [z*x*C-y*s, z*y*C+x*s, c+z*z*C]])

def load(name):
    p = PDBParser(QUIET=True)
    s = p.get_structure(name, os.path.join(PDBDIR, name+".pdb"))
    return s[0]  # first model

def his_coord_atoms(model):
    """Return dict resnum -> {'ND1':xyz,'NE2':xyz} for HIS residues."""
    d = {}
    ch = list(model.get_chains())[0]
    for res in ch:
        if res.get_resname() == "HIS":
            e = {}
            for an in ("ND1","NE2"):
                if an in res:
                    e[an] = res[an].get_coord().astype(float)
            if e:
                d[res.id[1]] = e
    return d

def pick_pair(hd, forced=None):
    """Pick coordinating His pair + specific N atoms. forced=(r1,r2) or closest."""
    if forced:
        r1, r2 = forced
        best = None
        for a1 in hd[r1]:
            for a2 in hd[r2]:
                dd = np.linalg.norm(hd[r1][a1]-hd[r2][a2])
                if best is None or dd < best[0]:
                    best = (dd, r1, a1, r2, a2)
        return best
    best = None
    for r1, r2 in itertools.combinations(sorted(hd), 2):
        for a1 in hd[r1]:
            for a2 in hd[r2]:
                dd = np.linalg.norm(hd[r1][a1]-hd[r2][a2])
                if best is None or dd < best[0]:
                    best = (dd, r1, a1, r2, a2)
    return best

def heavy_coords(model):
    xs = []
    for atom in model.get_atoms():
        el = atom.element.strip()
        if el == "H":
            continue
        xs.append(atom.get_coord().astype(float))
    return np.array(xs)

def build(name, forced=None, tag=""):
    model = load(name)
    hd = his_coord_atoms(model)
    dd, r1, a1n, r2, a2n = pick_pair(hd, forced)
    a1 = hd[r1][a1n]; a2 = hd[r2][a2n]
    ma = (a1+a2)/2.0
    u_mono = (a2-a1); u_mono/=np.linalg.norm(u_mono)
    # template occupied vertices t1=TET[0], t2=TET[1]
    t1, t2 = TET[0], TET[1]
    mt = (t1+t2)/2.0
    u_temp = (t2-t1); u_temp/=np.linalg.norm(u_temp)
    R_base = rot_align(u_mono, u_temp)
    strain = dd - VERT_SPACING  # induced-fit cost (N-N mismatch)

    heavy = heavy_coords(model)  # (N,3) monomer heavy atoms
    heavy_c = heavy - ma          # centered on N-midpoint

    def chainA(theta):
        Rj = rot_axis(u_temp, theta)
        Rf = Rj @ R_base
        A = (Rf @ heavy_c.T).T + mt
        return A, Rf
    def coordN(theta):
        Rj = rot_axis(u_temp, theta); Rf = Rj @ R_base
        n1 = Rf @ (a1-ma) + mt
        n2 = Rf @ (a2-ma) + mt
        return n1, n2

    # jigger scan
    angles = np.deg2rad(np.arange(0, 360, 5))
    results = []
    for th in angles:
        A,_ = chainA(th)
        B = (C2 @ A.T).T
        tree = cKDTree(B)
        pairs = tree.query_ball_point(A, 2.5)
        nclash = sum(len(p) for p in pairs)
        dmin = tree.query(A, k=1)[0].min()
        results.append((th, nclash, dmin))
    results = np.array(results)
    # minimize clash count, tiebreak by larger min-distance
    order = sorted(range(len(results)), key=lambda i:(results[i][1], -results[i][2]))
    best_i = order[0]
    best_th = results[best_i][0]
    best_clash = int(results[best_i][1])
    best_dmin = results[best_i][2]

    # realized geometry at best theta
    A, Rf = chainA(best_th)
    B = (C2 @ A.T).T
    n1, n2 = coordN(best_th)
    n3, n4 = C2@n1, C2@n2  # chain B coordinating N (C2 of chain A's)
    Ns = [n1, n2, n3, n4]
    dN = [np.linalg.norm(n) for n in Ns]  # Zn at origin
    angs = []
    for i,j in itertools.combinations(range(4),2):
        cA = np.dot(Ns[i],Ns[j])/(dN[i]*dN[j])
        angs.append(np.degrees(np.arccos(np.clip(cA,-1,1))))

    # write dimer PDB (no Zn), chain A + chain B, remap chain B to 'B'
    write_dimer(name, model, Rf, ma, mt, tag)

    return dict(name=name, tag=tag, r1=r1, a1n=a1n, r2=r2, a2n=a2n, dd=dd,
                strain=strain, best_th_deg=np.degrees(best_th), best_clash=best_clash,
                best_dmin=best_dmin, dN=dN, angs=angs,
                clash_curve=results)

def write_dimer(name, model, Rf, ma, mt, tag):
    """Apply transform to full model -> chain A; C2 -> chain B; write."""
    import copy
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model as BModel
    from Bio.PDB.Chain import Chain
    # transform chain A coords in place on a copy
    src = list(model.get_chains())[0]
    A = copy.deepcopy(src); A.id = "A"
    for atom in A.get_atoms():
        x = atom.get_coord().astype(float)
        atom.set_coord((Rf @ (x-ma)) + mt)
    B = copy.deepcopy(A); B.id = "B"
    for atom in B.get_atoms():
        atom.set_coord(C2 @ atom.get_coord().astype(float))
    st = Structure(name); m = BModel(0); st.add(m); m.add(A); m.add(B)
    io = PDBIO(); io.set_structure(st)
    out = os.path.join(OUT, f"{name}{tag}_dimer.pdb")
    io.save(out)
    return out

TARGETS = [
    ("UA_STRG.15519.1_1_193_357", (32,36), "_H32H36"),
    ("ES_STRG.14265.1_4_293_517", None, "_closest"),
    ("DK_STRG.14275.1_5_1279_1518", None, "_closest"),
    ("DK_STRG.12883.1_4_375_656", None, "_closest"),
]

if __name__ == "__main__":
    allres = []
    for nm, forced, tag in TARGETS:
        try:
            r = build(nm, forced, tag)
            allres.append(r)
            print(f"\n=== {nm}{tag} ===")
            print(f"  His pair: {r['r1']}({r['a1n']}) - {r['r2']}({r['a2n']}),  N-N = {r['dd']:.2f} A")
            print(f"  induced-fit strain (N-N minus {VERT_SPACING:.3f}) = {r['strain']:+.2f} A")
            print(f"  best jigger angle = {r['best_th_deg']:.0f} deg, min clash pairs(<2.5A) = {r['best_clash']}, min inter-chain dist = {r['best_dmin']:.2f} A")
            print(f"  4 N-Zn dist: " + ", ".join(f"{d:.2f}" for d in r['dN']))
            print(f"  6 N-Zn-N angles: " + ", ".join(f"{a:.0f}" for a in r['angs']))
            cc = r['clash_curve']
            print(f"  clash range over jigger: min={int(cc[:,1].min())} max={int(cc[:,1].max())}; #angles clash-free={int((cc[:,1]==0).sum())}/{len(cc)}")
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"FAIL {nm}: {e}")
    # save summary
    import json
    with open(os.path.join(OUT,"summary.txt"),"w") as f:
        for r in allres:
            f.write(f"{r['name']}{r['tag']}: pair {r['r1']}{r['a1n']}-{r['r2']}{r['a2n']} NN={r['dd']:.2f} strain={r['strain']:+.2f} "
                    f"bestang={r['best_th_deg']:.0f} clash={r['best_clash']} dmin={r['best_dmin']:.2f} "
                    f"dN={[round(x,2) for x in r['dN']]} angs={[round(x,0) for x in r['angs']]}\n")
    print("\nDONE")
