import pymupdf as fitz, json, numpy as np, csv, sys, collections
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pdf, raw, outd = sys.argv[1:4]; d = fitz.open(pdf); raw = json.load(open(raw))
def load(f):
    g = collections.OrderedDict()
    for r in csv.DictReader(open(f)): g.setdefault(int(r['id']), []).append(r)
    return g
cyc = load(f'{outd}/cyclic-loading.csv'); qs = load(f'{outd}/quasi-static.csv')
plan = [('4b', [(1, None)]), ('10a', [(2, None)])] + [(p, [(i, None) for i in ids]) for p, ids in
        [('5a', range(10, 13)), ('5b', range(13, 16)), ('6a', range(16, 19)), ('6b', range(19, 22))]]
fig, axs = plt.subplots(2, 3, figsize=(15, 9)); axs = axs.ravel()
cols = {0.384: 'tab:red', 0.128: 'tab:blue', 0.032: 'tab:green'}
for ax, (p, ids) in zip(axs, plan):
    m = raw[p]['_meta']; box = fitz.Rect(m['box'])
    pix = d[m['page']].get_pixmap(clip=box, dpi=400)
    img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)
    ax.imshow(img, extent=[0, m['xlim'], 0, m['ylim']], aspect='auto', alpha=0.8)
    for i, _ in ids:
        src = qs if p in ('4b', '10a') else cyc; rr = src[i]
        x = [float(r['stretch']) - 1 for r in rr]; y = [float(r['stress']) / 1e6 for r in rr]
        c = 'k' if p in ('4b', '10a') else cols[float(rr[0]['rate'])]
        ax.plot(x, y, 'o', ms=3, mfc='none', color=c, label=f"id {i}" + ('' if p in ('4b', '10a') else f", {rr[0]['rate']} /s"))
    ax.set_xlim(0, m['xlim']); ax.set_ylim(0, m['ylim']); ax.set_title(f'Fig. {p}: CSV points over the PDF figure')
    ax.set_xlabel('nominal strain = λ − 1'); ax.set_ylabel('nominal stress [MPa]'); ax.legend(fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig(f'{outd}/check_overlay.png', dpi=110)
