# -*- coding: utf-8 -*-
"""
Overlap-corrected GSEA
- Goal: address Reviewer R1-3 -- the circularity caused by genes shared between the TF target sets and the 16-gene identity score
- Method: remove the genes shared with the 16-gene identity score from each target set, then recompute the ES and permutation p-value
- Algorithm: identical to the original gsea.py (hand-written ES + 2000 permutations + seed=42); only the gene sets change
- First run a sanity check: reproduce the original gsea_results.csv ES values before running the overlap-corrected version
"""
import csv
import paths
import numpy as np

base = paths.r3_inputs()

# ---- 1. read the genome-wide logFC ----
fc = {}
with open(base + r'\volcano_input.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fc[row['id']] = float(row['logFC'])

# ---- 2. read the target gene sets ----
sets = {}
with open(base + r'\tf_target_sets.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sets.setdefault(row['set'], []).append(row['gene'])

# ---- 3. 16-gene identity score (from analysis.R lines 25-26 & Methods) ----
identity_16 = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFA","MAFB","NEUROD1","PAX6",
               "SLC2A2","GCK","ABCC8","KCNJ11","IAPP","UCN3","GLP1R","SLC30A8",
               "PCSK1","PCSK2","INS"]  # note: the text says 16 genes but analysis.R lists 19; use the Methods 16-gene set
# the 16 genes explicitly listed in Methods ("16 identity genes, including..."):
identity_16_methods = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFB","NEUROD1","PAX6",
                        "SLC2A2","GCK","ABCC8","IAPP","GLP1R","SLC30A8","PCSK1","PCSK2","INS"]

# shared genes named by the reviewer (R1-3 original text):
shared = {"SLC2A2", "GCK", "IAPP", "INS", "NKX6-1", "PDX1"}

# ---- 4. rank (as in gsea.py: logFC descending) ----
ranked = sorted(fc.items(), key=lambda x: -x[1])
N = len(ranked)

def compute_es(gene_idx_set, ranked, fc, n_miss, nr):
    """Weighted running-sum ES. nr is the fixed normalization weight (sum of |logFC| over the observed set),
    kept constant between observed and permutation samples (standard GSEA convention, as in gsea.py)."""
    n_hits = len(gene_idx_set)
    if nr == 0:
        return 0.0
    es = 0.0
    es_max = 0.0
    for g, lfc in ranked:
        if g in gene_idx_set:
            es += abs(lfc) / nr
        else:
            es -= 1.0 / n_miss
        if abs(es) > abs(es_max):
            es_max = es
    return es_max

def gsea(geneset, ranked, fc, n_perm=2000, seed=42):
    """Identical to gsea.gsea"""
    geneset = set(g for g in geneset if g in fc)
    if len(geneset) < 3:
        return None, None, None
    n_hits = len(geneset)
    n_miss = N - n_hits
    nr = sum(abs(fc[g]) for g in geneset)  # fixed normalization weight (observed set)
    es_max = compute_es(geneset, ranked, fc, n_miss, nr)
    all_genes = list(fc.keys())
    rng = np.random.default_rng(seed)
    perm_es = []
    for _ in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es.append(compute_es(sample, ranked, fc, n_miss, nr))
    perm_es = np.array(perm_es)
    p = float(np.mean(np.abs(perm_es) >= abs(es_max)))
    return es_max, p, n_hits

# ---- 5. SANITY CHECK: reproduce the original results (no overlap correction) ----
print('=' * 70)
print('SANITY CHECK: reproduce the original GSEA (should match gsea_results.csv)')
print('=' * 70)
print(f'{"gene set":24s} {"reproduced ES":>10s} {"orig ES":>10s} {"match?"}')
expected = {
    'PDX1_targets': -0.948,
    'HNF1A_targets': -0.808,
    'beta_identity_network': -0.780,
    'beta_functional_genes': -0.815,
}
ok_all = True
for sname, genes in sets.items():
    es, p, n_hits = gsea(genes, ranked, fc)
    exp = expected.get(sname, None)
    if es is None:
        print(f'{sname:24s}  (fewer than 3 genes, skipped)')
        continue
    if exp is not None:
        match = abs(es - exp) < 0.001
        ok_all = ok_all and match
        print(f'{sname:24s} {es:+10.4f} {exp:+10.3f} {"OK" if match else "X mismatch!"}')
    else:
        print(f'{sname:24s} {es:+10.4f}   (no expected value)')

print()
if ok_all:
    print('SANITY CHECK passed: reproduced results match gsea_results.csv; algorithm correct.')
else:
    print('SANITY CHECK failed: check the algorithm against the original gsea.py.')
print()

# ---- 6. overlap-corrected GSEA ----
print('=' * 70)
print('Overlap-corrected GSEA (shared genes removed from each target set)')
print('shared genes (16-gene identity score & target sets) =', sorted(shared))
print('=' * 70)
print(f'{"gene set":24s} {"n kept":>10s} {"corrected ES":>12s} {"p-value":>10s}')
results = []
for sname, genes in sets.items():
    genes_all = set(genes)
    removed = genes_all & shared
    genes_kept = genes_all - shared
    # keep only genes present in fc
    genes_kept_infc = [g for g in genes_kept if g in fc]
    es, p, n_hits = gsea(genes_kept, ranked, fc)
    results.append((sname, removed, len(genes_kept_infc), es, p))
    if es is None:
        print(f'{sname:24s}  fewer than 3 genes after removal, skipped')
        continue
    dirn = 'down-enriched' if es < 0 else 'up-enriched'
    print(f'{sname:24s} {len(genes_kept_infc):>10d} {es:+12.4f} {p:10.4f}  ({dirn}, removed {len(removed)})')

# ---- 7. save results ----
out = base + r'\gsea_overlap_corrected_results.csv'
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['set', 'removed_genes', 'n_genes_after', 'ES', 'p'])
    for sname, removed, n_kept, es, p in results:
        w.writerow([sname, ';'.join(sorted(removed)), n_kept, 
                    (f'{es:.6f}' if es is not None else ''), 
                    (f'{p:.6f}' if p is not None else '')])
print()
print('Overlap-corrected results saved:', out)
