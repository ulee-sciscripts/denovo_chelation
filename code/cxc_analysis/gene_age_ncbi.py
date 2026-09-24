#!/usr/bin/env python3
"""Gene age (phylostratum) per FBgn from NCBI Orthologs (gene/DATA/gene_orthologs.gz, downloaded
2026-09-23; D. mel rows only). Partner lineages from NCBI Taxonomy efetch. NCBI's arthropod ortholog
set has no non-arthropod partners, so the oldest stratum is 'Arthropoda or older'.
Writes gene_age_ncbi.json: {FBgn: stratum index 0..5}."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import gzip, json, re
import xml.etree.ElementTree as ET
D = _HERE
STRATA = ["D.mel only / no ortholog", "melanogaster subgroup/group", "Sophophora/Drosophila/Drosophilidae",
          "Diptera", "Insecta", "Arthropoda or older"]
def stratum(lin):
    L = lin.split("; ")
    if "melanogaster group" in L or "melanogaster subgroup" in L: return 1
    if "Drosophilidae" in L: return 2
    if "Diptera" in L: return 3
    if "Insecta" in L: return 4
    return 5
tax = {t.findtext("TaxId"): stratum((t.findtext("Lineage") or "") + "; " + (t.findtext("ScientificName") or ""))
       for t in ET.parse(f"{D}/data/ncbi_partner_taxonomy.xml").getroot().findall("Taxon")}
fb = {}
for line in gzip.open(f"{D}/data/ncbi_dmel_gene_info.gz", "rt"):
    if line.startswith("#"): continue
    f = line.split("\t"); m = re.search(r"FLYBASE:(FBgn\d+)", f[5])
    if m: fb[f[1]] = m.group(1)
age = {}
for line in gzip.open(f"{D}/data/ncbi_gene_orthologs_dmel.tsv.gz", "rt"):
    a, g1, _, b, g2 = line.rstrip("\n").split("\t")
    me, other = (g1, b) if a == "7227" else (g2, a)
    if me in fb and other in tax:
        age[fb[me]] = max(age.get(fb[me], 0), tax[other])
for g in fb.values(): age.setdefault(g, 0)
json.dump(age, open(f"{D}/data/gene_age_ncbi.json", "w"))
import collections; c = collections.Counter(age.values())
print(f"{len(age)} FBgn with NCBI GeneID;", {STRATA[k]: v for k, v in sorted(c.items())})
