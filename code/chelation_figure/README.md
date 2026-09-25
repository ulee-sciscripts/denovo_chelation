# chelation_figure — multimeric Zn²⁺ sequestration panel (aim2_chelation_panel)

Predicted homodimeric Zn-coordination sites for four de novo ORFs, each labelled with its
Metal3D dimer site-probability:

| panel | ORF | Metal3D (dimer) | model |
|---|---|---|---|
| A | ZMEG (UA_STRG.15519.1) | 0.95 | idealized symmetric C2 construct dimer (His32/His36) |
| B | DK14275 (DK_STRG.14275.1) | 0.85 | idealized symmetric C2 construct dimer (His54/His58) |
| C | ES14265 (ES_STRG.14265.1) | 0.79 | idealized symmetric C2 construct dimer (His40/His44) |
| D | O433 (FI_STRG.14643.1) | 0.99 | Boltz-2-derived His₂Cys₂ homodimer |

## 1. Rebuild the panel (self-contained, no PyMOL)
```
cd code/chelation_figure
python3 make_aim2_figure.py       # renders/*.png  ->  ../figures/aim2_chelation_panel.png (300 dpi)
```
`make_aim2_figure.py` (Pillow) crops each committed render in `renders/`, adds the panel
letter / ORF name / Metal3D score, and tiles the four. Bold font = Arial if present, else
Liberation/DejaVu Bold.

## 2. Regenerate the renders (needs open-source PyMOL)
The structures are committed under `structures/`, so the ray-traced renders are reproducible:
```
cd code/chelation_figure
pymol -cq render_dimer_panels.py      # structures/*.pdb -> renders/dimer_{zmeg,dk14275,es14265,fi14643}_downz.png
pymol -cq render_fi14643_his2cys2.py  # structures/FI_..._bestZn_model3.cif -> renders/dimer_fi14643_his2cys2_downz.png
```
Style: cartoon, coordinating His as red sticks, Zn sphere, white background, ray-traced 300 dpi.
(The panel uses the ZMEG/DK14275/ES14265 renders + the O433 His₂Cys₂ render.)

## 3. Structure provenance
- **Panels A–C (construct dimers).** `build_dimer.py` builds idealized symmetric-C2 homodimers:
  Zn at the origin with four coordinating imidazole N on tetrahedral vertices (Zn–N 2.05 Å),
  chain B the C2 image of chain A, rotationally jiggered to minimize inter-chain clash. The
  committed `structures/*_closest_dimer.pdb` / `*_H32H36_dimer.pdb` are its output. Metal3D was
  scored apo (Zn removed), site-restricted (boxes on the coordinating His Cα): ZMEG 0.946,
  DK14275 0.845, ES14265 0.786. `build_dimer.py` needs the monomer ESMFold models as input
  (external; set `MONOMER_PDBDIR` or place under `structures/monomers/`).
- **Panel D (O433).** `structures/FI_STRG_14643_1_8_2494_2775_bestZn_model3.cif` is the Boltz-2
  homodimer with an inter-chain His₂Cys₂ site (Cys62 + His89 per chain); Metal3D 0.9993.

## Files
- `make_aim2_figure.py` — assembles `../figures/aim2_chelation_panel.png` from `renders/`.
- `render_dimer_panels.py`, `render_fi14643_his2cys2.py` — PyMOL render scripts (need PyMOL).
- `build_dimer.py` — idealized-C2 homodimer construction (needs external monomer models).
- `structures/` — committed dimer PDBs + the O433 Boltz-2 CIF (render inputs).
- `renders/` — committed ray-traced PyMOL PNGs (panel inputs).
