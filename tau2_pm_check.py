# [S3.18] Tau-squared re-estimated by Paule-Mandel with Hartung-Knapp
# Supplementary 3 of the manuscript. See README.md for the full map.
"""[S3.18] Reestimativa de tau2 por Paule-Mandel com ajuste de Hartung-Knapp.
Testa se o tier consistente da descoberta e artefato do truncamento de DerSimonian-Laird.
Supplementary 3 of the manuscript. See README.md for the full map."""
import os, glob, argparse
import numpy as np, pandas as pd
from scipy import stats
import meta_pairs as mp

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.path.expanduser("~/syncom"))
ap.add_argument("--corte", default="2024-12-31")
a = ap.parse_args()
mm = os.path.join(a.root, "merged", "meta")

# --- estudos da descoberta, pelo mesmo arquivo de datas ---
cand = [f for f in glob.glob(os.path.join(a.root, "merged", "*.tsv"))
        if any(t in os.path.basename(f).lower() for t in ("date", "data"))]
if not cand:
    raise SystemExit("arquivo de datas nao encontrado em merged/: " +
                     str([os.path.basename(f) for f in glob.glob(os.path.join(a.root,'merged','*.tsv'))]))
dates = pd.read_csv(cand[0], sep="\t", dtype=str)
print("datas:", os.path.basename(cand[0]), "| colunas:", list(dates.columns))
cs = [c for c in dates.columns if "stud" in c.lower() or "acc" in c.lower()][0]
cd = [c for c in dates.columns if c != cs][0]
dt = dict(zip(dates[cs], dates[cd]))
desc = {s for s, v in dt.items() if isinstance(v, str) and v <= a.corte}
print(f"estudos na descoberta: {len(desc)}")

d, _qc = mp.load(a.root, False)
if "study" not in d.columns:
    raise SystemExit("load() sem coluna study; colunas: " + str(list(d.columns)))
d = d[d.study.isin(desc)]
d = d[np.isfinite(d.z) & (d.se > 0)].copy()
d["pair"] = d.genus_a + "|" + d.genus_b

D0 = pd.read_csv(os.path.join(mm, "validacao_temporal", "descoberta.tsv.gz"), sep="\t")
min_k = int(D0.k.min()); print(f"min_k da descoberta: {min_k}")
kk = d.groupby("pair").z.transform("size")
d = d[kk >= min_k].copy()

code, pares = pd.factorize(d.pair.values)
z = d.z.values.astype(float); v = (d.se.values.astype(float))**2
k = np.bincount(code).astype(float); npair = len(pares)
print(f"pares: {npair} | observacoes: {len(z)}")

def somas(tau2):
    w = 1.0/(v + tau2[code])
    sw = np.bincount(code, w); swz = np.bincount(code, w*z)
    mu = swz/sw
    q = np.bincount(code, w*(z - mu[code])**2)
    return w, sw, mu, q

# Paule-Mandel por bisseccao vetorizada: resolve Q(tau2) = k-1
lo = np.zeros(npair); hi = np.full(npair, 5.0)
for _ in range(40):
    _, _, _, q = somas(hi)
    alto = q > (k - 1)
    if not alto.any(): break
    hi[alto] *= 2
for _ in range(60):
    mid = 0.5*(lo + hi)
    _, _, _, q = somas(mid)
    acima = q > (k - 1)
    lo = np.where(acima, mid, lo); hi = np.where(acima, hi, mid)
tau2_pm = 0.5*(lo + hi)

w, sw, mu, q = somas(tau2_pm)
se_hk = np.sqrt(np.maximum(q/((k - 1)*sw), 1e-12))          # Hartung-Knapp
df = np.maximum(k - 1, 1)
p_pm = 2*stats.t.sf(np.abs(mu/se_hk), df)
half = stats.t.ppf(0.975, df)*np.sqrt(tau2_pm + se_hk**2)
pi_lo, pi_hi = mu - half, mu + half

o = np.argsort(p_pm); m = npair; fdr = np.empty(m)
fdr[o] = np.minimum.accumulate((p_pm[o]*m/np.arange(1, m+1))[::-1])[::-1]
fdr = np.minimum(fdr, 1)

sig = fdr < 0.05
cls = np.where(~sig, "ns", np.where(mu > 0, "associado_pos", "associado_neg"))
cls = np.where(sig & (pi_lo > 0), "consistente_pos", cls)
cls = np.where(sig & (pi_hi < 0), "consistente_neg", cls)

P = pd.DataFrame({"pair": pares, "k": k.astype(int), "mu_pm": mu, "se_hk": se_hk,
                  "tau2_pm": tau2_pm, "pi_lo_pm": pi_lo, "pi_hi_pm": pi_hi,
                  "fdr_pm": fdr, "classe_pm": cls})
P.to_csv(os.path.join(mm, "validacao_temporal", "descoberta_pm.tsv.gz"),
         sep="\t", index=False, float_format="%.5g")

print("\n=== tau2: DL contra PM (pares em comum) ===")
cmp = P.merge(D0.assign(pair=D0.genus_a + "|" + D0.genus_b)[["pair","tau2","classe","fdr"]], on="pair")
print(f"tau2 == 0: DL {100*(cmp.tau2<=1e-8).mean():.1f}%  |  PM {100*(cmp.tau2_pm<=1e-8).mean():.1f}%")
print(f"tau2 mediano: DL {cmp.tau2.median():.4f}  |  PM {cmp.tau2_pm.median():.4f}")
print("\nclasse DL -> PM:")
print(pd.crosstab(cmp.classe, cmp.classe_pm).to_string())

V = pd.read_csv(os.path.join(mm, "validacao_temporal", "validacao.tsv"), sep="\t")
J = V.merge(P[["pair","classe_pm","tau2_pm","fdr_pm"]], on="pair", how="inner")
J["grupo_pm"] = np.where(J.classe_pm.str.startswith("consistente"), "consistente",
                  np.where(J.classe_pm.str.startswith("associado"), "associado", "ns"))
print("\n=== replicacao por classificacao ===")
for nome, col in (("DL (atual)", "grupo"), ("PM + HK", "grupo_pm")):
    print(f"\n{nome}:")
    for g, s in J[J[col] != "ns"].groupby(col):
        print(f"  {g:<14} n={len(s):>5}  replica={100*(s.replicado==True).mean():.0f}%  "
              f"k_val mediano={s.k_val.median():.0f}  |r|val={s.r_mu_val.abs().median():.3f}")
