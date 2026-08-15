"""
Altschul-Erikson (1985) dinucleotide-preserving sequence shuffle via a random
Eulerian path on the doublet multigraph. Preserves EXACT dinucleotide counts and
the first/last symbol. Verified by direct dinucleotide-count comparison.

Reference: Altschul SF, Erikson BW (1985) Mol Biol Evol 2:526-538;
Eulerian-path formulation (Kandel et al. 1996). Alphabet here is A/C/G/T.
"""
import numpy as np


def _dinuc_counts(seq):
    d = {}
    for i in range(len(seq) - 1):
        k = seq[i:i + 2]
        d[k] = d.get(k, 0) + 1
    return d


def altschul_erikson_shuffle(seq, rng, max_tries=100):
    """Return a shuffle of `seq` with identical dinucleotide composition and the
    same first & last symbol. `rng` is a numpy Generator. Falls back to the input
    only for degenerate seqs (len<3 or a vertex with no exit that isn't last)."""
    n = len(seq)
    if n < 3:
        return seq
    symbols = sorted(set(seq))
    last = seq[-1]
    first = seq[0]
    # outgoing edges (as target symbols) per vertex, in original order
    edges = {s: [] for s in symbols}
    for i in range(n - 1):
        edges[seq[i]].append(seq[i + 1])

    for _ in range(max_tries):
        # 1) choose, for each vertex != last, one outgoing edge to be traversed LAST
        last_choice = {}
        feasible = True
        for v in symbols:
            if v == last:
                continue
            if not edges[v]:
                feasible = False
                break
            last_choice[v] = int(rng.integers(len(edges[v])))
        if not feasible:
            return seq  # degenerate; caller may skip
        # 2) the chosen "last edges" must form an arborescence rooted at `last`
        good = True
        for v in symbols:
            if v == last:
                continue
            seen = set()
            cur = v
            while cur != last:
                if cur in seen or cur not in last_choice:
                    good = False
                    break
                seen.add(cur)
                cur = edges[cur][last_choice[cur]]
            if not good:
                break
        if not good:
            continue
        # 3) per vertex, random-permute the non-last edges, append the last edge
        new_edges = {}
        for v in symbols:
            elist = edges[v][:]
            if v == last or v not in last_choice:
                rng.shuffle(elist)
                new_edges[v] = elist
            else:
                li = last_choice[v]
                le = elist[li]
                rest = elist[:li] + elist[li + 1:]
                rng.shuffle(rest)
                new_edges[v] = rest + [le]
        # 4) walk the Eulerian path from `first`
        out = [first]
        ptr = {v: 0 for v in symbols}
        cur = first
        ok = True
        for _ in range(n - 1):
            lst = new_edges[cur]
            p = ptr[cur]
            if p >= len(lst):
                ok = False
                break
            nxt = lst[p]
            ptr[cur] = p + 1
            out.append(nxt)
            cur = nxt
        if ok and len(out) == n:
            return "".join(out)
    return seq  # give up (extremely rare); caller may skip


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    bad = 0
    for _ in range(2000):
        L = int(rng.integers(30, 200))
        s = "".join(rng.choice(list("ACGT"), size=L))
        sh = altschul_erikson_shuffle(s, rng)
        if _dinuc_counts(s) != _dinuc_counts(sh):
            bad += 1
        if s[0] != sh[0] or s[-1] != sh[-1] or len(s) != len(sh):
            bad += 1
    print("dinuc-count / endpoint mismatches over 2000 random seqs:", bad)
