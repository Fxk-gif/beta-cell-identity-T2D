import csv, os
import paths
import numpy as np
from scipy import stats

base = paths.r3_inputs()

# read the genome-wide logFC
fc = {}
with open(base + r'\volcano_input.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fc[row['id']] = float(row['logFC'])

# read the target gene sets
sets = {}
with open(base + r'\tf_target_sets.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sets.setdefault(row['set'], []).append(row['gene'])

# sort by logFC descending (standard GSEA: positive logFC first)
ranked = sorted(fc.items(), key=lambda x: -x[1])  # (gene, logFC) descending
N = len(ranked)

def compute_es(gene_idx_set, ranked, fc, n_miss, nr):
    """Weighted running-sum ES. nr is the fixed normalization weight (sum of |logFC| over the observed set),
    kept constant between observed and permutation samples (standard GSEA convention)."""
    n_hits = len(gene_idx_set)
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

def gsea(geneset, ranked, n_perm=2000):
    """Hand-written GSEA; returns the ES and p-value"""
    geneset = set(g for g in geneset if g in fc)
    if len(geneset) < 3:
        return None, None, None
    n_hits = len(geneset)
    n_miss = N - n_hits
    nr = sum(abs(fc[g]) for g in geneset)  # fixed normalization weight (observed set)
    # observed ES
    es_max = compute_es(geneset, ranked, fc, n_miss, nr)
    # permutation test
    all_genes = list(fc.keys())
    perm_es = []
    rng = np.random.default_rng(42)
    for _ in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es.append(compute_es(sample, ranked, fc, n_miss, nr))
    perm_es = np.array(perm_es)
    p = float(np.mean(np.abs(perm_es) >= abs(es_max)))
    if es_max >= 0:
        p_up = float(np.mean(perm_es >= es_max))
    else:
        p_up = float(np.mean(perm_es <= es_max))
    return es_max, p, n_hits

print('=== GSEA enrichment test (genome-wide logFC ranking) ===')
print(f'{"gene set":22s} {"ES":>8s} {"p(two-sided)":>10s} {"n_hits":>7s}  direction')
results = []
for sname, genes in sets.items():
    es, p, n_hits = gsea(genes, ranked)
    if es is None:
        continue
    direction = 'down-enriched' if es < 0 else 'up-enriched'
    results.append((sname, es, p, n_hits))
    print(f'{sname:22s} {es:+8.3f} {p:10.4f} {n_hits:7d}  {direction}')

# save results
with open(base + r'\gsea_results.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['set', 'ES', 'p', 'n_hits'])
    for r in results:
        w.writerow(list(r))
print('\nGSEA results saved: gsea_results.csv')
