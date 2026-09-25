#!/usr/bin/env python3
"""Build the Aim 2 four-panel figure for the research statement from the dimer
renders. Each panel: alpha-cropped structure render on white, bold panel letter
(top-left), monomer short-name (top-right), Metal3D score (bottom-right).
Style matches the statement's other figures (Arial-bold labels, white ground).

Panels (score = Metal3D dimer site-probability):
  A ZMEG    0.95  (site-restricted construct dimer; RESULTS.md 0.946)
  B DK14275 0.85  (site-restricted construct dimer; 0.845)
  C ES14265 0.79  (site-restricted construct dimer; 0.786)
  D FI14643 0.99  (Boltz-derived His2Cys2 homodimer; metal3d_dimer_consistent 0.9993)

Run: python3 make_aim2_figure.py   ->  renders/aim2_chelation_panel.png
"""
import os
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
RND = os.path.join(HERE, "renders")
FIGS = os.path.join(HERE, "..", "figures")
# Bold sans-serif label font; prefer Arial (original render), fall back to bundled equivalents.
ARIAL_BD = next((p for p in [
    "/home/unjin/.local/share/fonts/arialbd.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
] if os.path.exists(p)), "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")

# (file, letter, name, score-string)
# (file, letter, name, score, zoom)  zoom>1 = crop tighter around the dense core
PANELS = [
    ("dimer_zmeg_downz.png",             "A", "ZMEG",    "0.95", 1.0),
    ("dimer_dk14275_downz.png",          "B", "DK14275", "0.85", 1.45),  # thin diagonal -> zoom in
    ("dimer_es14265_downz.png",          "C", "ES14265", "0.79", 1.0),
    ("dimer_fi14643_his2cys2_downz.png", "D", "O433", "0.99", 1.0),   # Boltz His2Cys2 model
]

PAD = 24            # px white margin around cropped content inside each cell
CELL_W, CELL_H = 1500, 1150   # cell aspect ~ content (landscape helices)
GAP = 20            # gap between cells
MARGIN = 12         # outer margin

def content_crop(im, pad=18):
    """Crop to the structure's bounding box. The renders have a WHITE (opaque)
    background, so detect non-white RGB pixels (not alpha). A small pad keeps the
    structure off the panel edge."""
    import numpy as np
    im = im.convert("RGB")
    arr = np.asarray(im)
    nonwhite = (arr < 245).any(axis=-1)            # any channel below near-white
    if not nonwhite.any():
        return im
    ys, xs = np.where(nonwhite)
    l, r = max(0, xs.min() - pad), min(im.width,  xs.max() + 1 + pad)
    t, b = max(0, ys.min() - pad), min(im.height, ys.max() + 1 + pad)
    return im.crop((l, t, r, b))

def zoom_core(im, zoom):
    """Crop a centered box (about the pixel centroid = the dense core) covering
    1/zoom of each dimension, so fit_into upscales it -> structure appears larger.
    Trims the far tail tips of elongated panels (e.g. the DK14275 diagonal)."""
    import numpy as np
    if zoom <= 1.0:
        return im
    arr = np.asarray(im.convert("RGB"))
    nonwhite = (arr < 245).any(axis=-1)
    ys, xs = np.where(nonwhite)
    cy, cx = ys.mean(), xs.mean()                 # centroid weights the dense region
    hw, hh = im.width / (2 * zoom), im.height / (2 * zoom)
    l = int(max(0, min(cx - hw, im.width  - 2 * hw)))
    t = int(max(0, min(cy - hh, im.height - 2 * hh)))
    return im.crop((l, t, int(l + 2 * hw), int(t + 2 * hh)))

def fit_into(content, cw, ch, pad):
    """White cell (cw x ch); paste aspect-preserved content centered, pad margin."""
    cell = Image.new("RGBA", (cw, ch), (255, 255, 255, 255))
    aw, ah = cw - 2 * pad, ch - 2 * pad
    cr = content.width / content.height
    ar = aw / ah
    if cr > ar:
        nw, nh = aw, max(1, round(aw / cr))
    else:
        nh, nw = ah, max(1, round(ah * cr))
    content = content.resize((nw, nh), Image.LANCZOS)
    ox, oy = pad + (aw - nw) // 2, pad + (ah - nh) // 2
    cell.alpha_composite(content, (ox, oy))
    return cell

def draw_labels(cell, letter, name, score):
    d = ImageDraw.Draw(cell)
    W, H = cell.size
    f_letter = ImageFont.truetype(ARIAL_BD, 150)
    f_name   = ImageFont.truetype(ARIAL_BD, 92)
    f_slbl   = ImageFont.truetype(ARIAL_BD, 52)
    f_score  = ImageFont.truetype(ARIAL_BD, 110)
    m = 26
    # panel letter, top-left
    d.text((m, m), letter, font=f_letter, fill=(0, 0, 0, 255))
    # monomer name, top-right
    tb = d.textbbox((0, 0), name, font=f_name)
    d.text((W - m - (tb[2] - tb[0]), m + 14), name, font=f_name, fill=(0, 0, 0, 255))
    # Metal3D score, bottom-right (small caps label above bold number)
    sb = d.textbbox((0, 0), score, font=f_score)
    lb = d.textbbox((0, 0), "Metal3D", font=f_slbl)
    sx = W - m - (sb[2] - sb[0])
    sy = H - m - (sb[3] - sb[1]) - 8
    lx = W - m - (lb[2] - lb[0])
    ly = sy - (lb[3] - lb[1]) - 16
    d.text((lx, ly), "Metal3D", font=f_slbl, fill=(90, 90, 90, 255))
    d.text((sx, sy), score, font=f_score, fill=(0, 0, 0, 255))
    return cell

def main():
    cells = []
    for fn, letter, name, score, zoom in PANELS:
        im = Image.open(os.path.join(RND, fn))
        content = zoom_core(content_crop(im), zoom).convert("RGBA")
        cell = fit_into(content, CELL_W, CELL_H, PAD)
        cell = draw_labels(cell, letter, name, score)
        cells.append(cell)

    W = MARGIN * 2 + CELL_W * 2 + GAP
    H = MARGIN * 2 + CELL_H * 2 + GAP
    canvas = Image.new("RGBA", (W, H), (255, 255, 255, 255))
    pos = [(MARGIN, MARGIN),
           (MARGIN + CELL_W + GAP, MARGIN),
           (MARGIN, MARGIN + CELL_H + GAP),
           (MARGIN + CELL_W + GAP, MARGIN + CELL_H + GAP)]
    for cell, p in zip(cells, pos):
        canvas.alpha_composite(cell, p)

    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, "aim2_chelation_panel.png")
    canvas.convert("RGB").save(out, dpi=(300, 300))
    print("wrote", out, canvas.size)

if __name__ == "__main__":
    main()
