import csv, os, sys
from collections import defaultdict

BASE = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\geo_data\processed'
ANN = os.path.join(BASE, 'gpl570_probe2symbol.csv')

# 读取注释: probe_id -> gene_symbol (取 SYMBOL 非空)
probe2sym = {}
with open(ANN, encoding='utf-8') as f:
    r = csv.DictReader(f)
    cols = r.fieldnames
    # hgu133plus2SYMBOL toTable: probe_id, symbol
    pcol = 'probe_id' if 'probe_id' in cols else cols[0]
    scol = 'symbol' if 'symbol' in cols else cols[-1]
    for row in r:
        pid = (row.get(pcol) or '').strip()
        sym = (row.get(scol) or '').strip()
        if pid and sym and sym not in ('', 'NA', '---'):
            probe2sym[pid] = sym
print('annotated probes:', len(probe2sym))

def collapse(gse):
    src = os.path.join(BASE, f'{gse}_expr_probe.csv')
    out_expr = os.path.join(BASE, f'{gse}_expr_gene.csv')
    if not os.path.exists(src):
        print(gse, 'missing', src); return
    with open(src, encoding='utf-8') as f:
        rd = csv.reader(f)
        header = next(rd)
        samples = header[1:]
        best = {}  # gene -> (row, mean_expr)
        n_total = n_mapped = 0
        for row in rd:
            n_total += 1
            pid = row[0]
            sym = probe2sym.get(pid)
            if not sym:
                continue
            n_mapped += 1
            vals = [float(x) for x in row[1:]]
            mean_v = sum(vals) / len(vals)
            if sym not in best or mean_v > best[sym][1]:
                best[sym] = (row, mean_v)
    genes = sorted(best.keys())
    with open(out_expr, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['Gene.Symbol'] + samples)
        for g in genes:
            w.writerow([g] + best[g][0][1:])
    print(f'{gse}: {n_total} probes -> {n_mapped} mapped -> {len(genes)} unique genes -> {out_expr}')

for gse in ['GSE76894', 'GSE76895']:
    collapse(gse)
print('DONE')
