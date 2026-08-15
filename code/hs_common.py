"""
Shared helpers for the de novo His-spacing / chelation analysis.
Pure-Python + numpy. No plotting here (plotting is base R).
"""
import os
import gzip
import numpy as np

# Raw-input root, resolved relative to this file (code/ -> ../data/raw).
# fetch_data.sh downloads the heavy raw inputs (Grandchamp neORF FASTAs, FlyBase
# proteome + GTF, UniProt metal-binders, UCSC dm6 genome + RepeatMasker) here.
SCRATCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "data", "raw")

# ---- standard genetic code -------------------------------------------------
_BASES = "TCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON2AA = {}
_i = 0
for a in _BASES:
    for b in _BASES:
        for c in _BASES:
            CODON2AA[a + b + c] = _AA[_i]
            _i += 1

AAS = list("ACDEFGHIKLMNPQRSTVWY")  # 20 standard aa (order for composition tables)


def read_fasta(path):
    """Yield (header, seq) from a plain or .gz FASTA."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        name, chunks = None, []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def translate(dna, truncate_at_stop=True):
    """Translate frame 0. If truncate_at_stop, stop at first '*' and drop it.
    Returns (peptide, coding_dna) where coding_dna is the nt used (multiple of 3,
    excluding the terminating stop codon when truncated)."""
    n = (len(dna) // 3) * 3
    aa = []
    for i in range(0, n, 3):
        codon = dna[i:i + 3]
        res = CODON2AA.get(codon, "X")
        if res == "*" and truncate_at_stop:
            return "".join(aa), dna[:i]
        aa.append(res)
    return "".join(aa), dna[:n]


def translate_full(dna):
    """Translate frame 0 over the whole sequence; stop codons become '*'.
    Length = len(dna)//3. Used for the dinucleotide-null (positions kept fixed)."""
    n = (len(dna) // 3) * 3
    out = []
    for i in range(0, n, 3):
        out.append(CODON2AA.get(dna[i:i + 3], "X"))
    return "".join(out)


def load_denovo(paths, min_aa=10):
    """Translate the seven de novo nt FASTAs, frame 0, truncate at first stop,
    keep peptides >= min_aa, dedupe on exact protein string (first kept).
    Returns list of dicts: {id, prot, dna} where dna is the coding nt (peptide*3)."""
    seen = set()
    out = []
    for p in paths:
        for hdr, seq in read_fasta(p):
            seq = seq.upper()
            prot, cds = translate(seq, truncate_at_stop=True)
            if len(prot) < min_aa:
                continue
            if prot in seen:
                continue
            seen.add(prot)
            out.append({"id": hdr.split()[0], "prot": prot, "dna": cds})
    return out


def load_proteins(path, min_aa=10, drop_interior_stop=True):
    """Load a protein FASTA (optionally .gz); keep >= min_aa; drop entries whose
    interior contains '*'. Returns list of peptide strings."""
    out = []
    for hdr, seq in read_fasta(path):
        seq = seq.strip().upper()
        if seq.endswith("*"):
            seq = seq[:-1]
        if len(seq) < min_aa:
            continue
        if drop_interior_stop and "*" in seq:
            continue
        out.append(seq)
    return out


# ---- His-spacing counting --------------------------------------------------
def residue_positions(prot, res="H"):
    return np.frombuffer(prot.encode(), dtype=np.uint8) == ord(res)


def count_spacings(prot, dmax=6, res="H"):
    """Return array c[1..dmax]: number of ordered start positions i with
    prot[i]==res and prot[i+d]==res, for d=1..dmax."""
    mask = residue_positions(prot, res)
    L = mask.size
    c = np.zeros(dmax + 1, dtype=np.int64)  # index by d
    for d in range(1, dmax + 1):
        if L > d:
            c[d] = int(np.count_nonzero(mask[:-d] & mask[d:]))
    return c[1:]  # d=1..dmax


def observed_spacing_totals(prots, dmax=6, res="H"):
    tot = np.zeros(dmax, dtype=np.int64)
    for p in prots:
        tot += count_spacings(p, dmax, res)
    return tot


def analytic_expected_totals(prots, dmax=6, res="H"):
    """Exact scramble expectation E_d = sum (L-d) h(h-1)/(L(L-1))."""
    E = np.zeros(dmax, dtype=np.float64)
    for p in prots:
        L = len(p)
        h = p.count(res)
        if L < 2 or h < 2:
            continue
        base = h * (h - 1) / (L * (L - 1))
        for d in range(1, dmax + 1):
            if L > d:
                E[d - 1] += (L - d) * base
    return E


def mc_expected_totals(prots, dmax=6, res="H", nrep=50, seed=0):
    """Monte-Carlo scramble expectation: for each protein place h residues at
    uniformly random positions among L, count spacings, average over nrep.
    Vectorised per (protein) using boolean masks."""
    rng = np.random.default_rng(seed)
    acc = np.zeros(dmax, dtype=np.float64)
    for p in prots:
        L = len(p)
        h = p.count(res)
        if L < 2 or h < 2:
            continue
        for _ in range(nrep):
            pos = rng.choice(L, size=h, replace=False)
            m = np.zeros(L, dtype=bool)
            m[pos] = True
            for d in range(1, dmax + 1):
                if L > d:
                    acc[d - 1] += np.count_nonzero(m[:-d] & m[d:])
    return acc / nrep


def fold_profile(prots, dmax=6, res="H", method="analytic", nrep=50, seed=0):
    obs = observed_spacing_totals(prots, dmax, res).astype(float)
    if method == "analytic":
        exp = analytic_expected_totals(prots, dmax, res)
    else:
        exp = mc_expected_totals(prots, dmax, res, nrep=nrep, seed=seed)
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = np.where(exp > 0, obs / exp, np.nan)
    return fold, obs, exp
