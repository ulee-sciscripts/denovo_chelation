# chelation_figure — multimeric Zn²⁺ sequestration panel (aim2_chelation_panel)

Assembles the four-panel figure of predicted homodimeric Zn-coordination sites for four
de novo ORFs, each labelled with its Metal3D dimer site-probability:

| panel | ORF | Metal3D (dimer) | model |
|---|---|---|---|
| A | ZMEG (UA_STRG.15519.1) | 0.95 | idealized symmetric C2 construct dimer (His32/His36) |
| B | DK14275 (DK_STRG.14275.1) | 0.85 | idealized symmetric C2 construct dimer (His54/His58) |
| C | ES14265 (ES_STRG.14265.1) | 0.79 | idealized symmetric C2 construct dimer (His40/His44) |
| D | O433 (FI_STRG.14643.1) | 0.99 | Boltz-2-derived His₂Cys₂ homodimer |

## Reproducing the figure (self-contained)
```
cd code/chelation_figure
python3 make_aim2_figure.py     # -> ../figures/aim2_chelation_panel.png (300 dpi)
```
`make_aim2_figure.py` (Pillow) crops each committed render in `renders/` to its structure,
adds the panel letter / ORF name / Metal3D score, and tiles the four into one panel. The four
input renders are committed here, so the figure rebuilds without PyMOL. A bold sans-serif font
is used for labels (Arial if present, else Liberation/DejaVu Bold).

## Upstream render pipeline (NOT committed — needs PyMOL + structures)
The `renders/dimer_*_downz.png` inputs are ray-traced PyMOL images (cartoon, coordinating
His shown as red sticks, Zn sphere; white background, 300 dpi). They are produced from
homodimer structures by `render_dimer_panels.py` / `render_fi14643_his2cys2.py` (open-source
PyMOL), which live with the structural pipeline in the project's chelation-biochem working set:
- Panels A–C: idealized symmetric C2 homodimers — Zn placed at the origin with four N on
  tetrahedral vertices (2.05 Å), chain B the C2 image of chain A, rotationally jiggered to
  minimize inter-chain clash. Metal3D scored site-restricted (boxes on the coordinating His
  Cα only), apo (Zn removed): ZMEG 0.946, DK14275 0.845, ES14265 0.786.
- Panel D: Boltz-2 homodimer of O433/FI14643 with an inter-chain His₂Cys₂ site (Cys62 + His89
  per chain); Metal3D 0.9993.

Only the assembler + its committed input renders are included here; regenerating the renders
themselves requires the external structures and PyMOL.

## Files
- `make_aim2_figure.py` — assembles `../figures/aim2_chelation_panel.png` from `renders/`.
- `renders/dimer_{zmeg,dk14275,es14265,fi14643_his2cys2}_downz.png` — committed PyMOL renders (inputs).
