#!/usr/bin/env python3
# [S1.12] Sample curation
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Aplica a curadoria (curadoria_manual.tsv) e varre todos os estudos atras de amostras fora da rizosfera.

curadoria_manual.tsv: study, tipo (fator|excluir), chave (atributo do BioSample), metodo, arg, nome
  metodo attr : o proprio valor do atributo
         tok  : parte n (arg, a partir de 1) do valor dividido por _ - . espaco
         re   : primeiro grupo da expressao regular (arg); exclusao = amostras em que a regex casa
Varredura automatica (todos os estudos incluidos): uma chave cujos valores variam entre as amostras
  e em que um padrao de outro compartimento (bulk, unplanted, grain, seed, negative control, blank,
  endosphere, leaf, shoot, nodule, tokens BS junto com RS) casa com parte (< 90%) das amostras
  -> amostras marcadas para exclusao (origem "auto"; conferir em excluir_amostras.tsv)

Saidas (<root>/merged/curadoria_final/):
  attrs_curado/        atributos + fatores derivados (cur_<nome>)  -> effect_sizes.py --attrs
  design_curado.tsv    study, fatores                               -> effect_sizes.py --design
  excluir_amostras.tsv study, run, origem, motivo                   -> effect_sizes.py --excluir
  resumo.txt
Uso: LD_PRELOAD=... python3 apply_curation.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs \
        --curadoria curadoria_manual.tsv
"""
import argparse, csv, gzip, os, re, shutil, sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import load_attrs

OTHER = re.compile(r"\bbulk\b|bulk[_ ]?soil|unplanted|non[- ]?rhizospher|grains?[ _](tissue|sample)|\bseeds?\b|(extraction|pcr|kit|reagent|sequencing)[ _](negative[ _])?(control|blank)|\bNTC\b|"
                   r"\bblank\b|mock[ _]community|zymobiomics|endospher|root[ _]interior|\bleaf\b|\bleaves\b|"
                   r"phyllospher|\bnodules?\b|stem[ _](tissue|sample)|shoot[ _](tissue|sample)", re.I)
SPLIT = re.compile(r"[_\-\.\s]+")


def samples_of(md, s, run2bs):
    with gzip.open(os.path.join(md, "genus", f"{s}.tsv.gz"), "rt") as fh:
        fh.readline()
        runs = [l.split("\t", 1)[0] for l in fh]
    return runs, [run2bs.get(r, r) for r in runs]


def derive(val, metodo, arg):
    if val is None:
        return None
    if metodo == "attr":
        return val
    if metodo == "tok":
        t = [x for x in SPLIT.split(val.strip()) if x]
        p = int(arg) - 1
        return t[p] if 0 <= p < len(t) else None
    m = re.search(arg, val)
    return (m.group(1) if m.groups() else m.group(0)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--curadoria", default="curadoria_manual.tsv")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); out = os.path.join(md, "curadoria_final")
    ad = os.path.join(out, "attrs_curado"); os.makedirs(ad, exist_ok=True)
    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    inc = [r["study"] for r in csv.DictReader(open(os.path.join(md, "study_qc.tsv")), delimiter="\t")
           if r["status"] == "incluido"]
    spec = defaultdict(list)
    for r in csv.DictReader(open(a.curadoria), delimiter="\t"):
        spec[r["study"]].append(r)
    for f in os.listdir(a.attrs):
        shutil.copy(os.path.join(a.attrs, f), os.path.join(ad, f))

    excl, design, L = [], [], []
    for s in sorted(inc):
        runs, bss = samples_of(md, s, run2bs)
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        n = len(runs)
        # varredura automatica
        keys = Counter(k for b in bss for k in attrs.get(b, {}))
        auto = set()
        for k in keys:
            vals = [attrs.get(b, {}).get(k, "") for b in bss]
            if len(set(vals)) < 2:
                continue
            hit = [bool(OTHER.search(v)) for v in vals]
            toks = [set(SPLIT.split(v)) for v in vals]
            bs = [("BS" in t) for t in toks]
            if any("RS" in t for t in toks) and 0 < sum(bs) < n:
                hit = [h or b for h, b in zip(hit, bs)]
            if 0 < sum(hit) < 0.9 * n:
                for r, b, h, v in zip(runs, bss, hit, vals):
                    if h:
                        m = OTHER.search(v)
                        ctx = v[max(0, m.start() - 30):m.end() + 30] if m else v[:60]
                        auto.add(r); excl.append((s, r, "auto", f"{k}: ...{ctx}..."))
        # curadoria manual
        new, fat, nman = [], [], 0
        for sp in spec.get(s, []):
            for r, b in zip(runs, bss):
                v = attrs.get(b, {}).get(sp["chave"])
                if sp["tipo"] == "excluir":
                    if v is not None and re.search(sp["arg"], v) and r not in auto:
                        excl.append((s, r, "manual", sp.get("nota", "")[:80])); nman += 1
                else:
                    d = derive(v, sp["metodo"], sp["arg"])
                    if d is not None:
                        new.append((b, f"cur_{sp['nome']}", d))
            if sp["tipo"] == "fator":
                fat.append(f"cur_{sp['nome']}")
        if new:
            p = os.path.join(ad, f"{s}.tsv")
            hdr = not os.path.exists(p) or os.path.getsize(p) == 0
            with open(p, "a") as fo:
                if hdr:
                    fo.write("sample\tkey\tvalue\n")
                for b, k, v in sorted(set(new)):
                    fo.write(f"{b}\t{k}\t{v}\n")
        if fat:
            design.append((s, ";".join(dict.fromkeys(fat))))
        ne = len({r for (st, r, _, _) in excl if st == s})
        if ne or fat:
            L.append(f"{s}: n={n} | excluidas={ne} (auto {len(auto)}, manual {nman}) | fatores: {', '.join(dict.fromkeys(fat)) or '-'}")

    with open(os.path.join(out, "design_curado.tsv"), "w") as fo:
        fo.write("study\tfatores\n")
        for s, f in design:
            fo.write(f"{s}\t{f}\n")
    seen = set()
    with open(os.path.join(out, "excluir_amostras.tsv"), "w") as fo:
        fo.write("study\trun\torigem\tmotivo\n")
        for s, r, o, m in excl:
            if (s, r) not in seen:
                seen.add((s, r)); fo.write(f"{s}\t{r}\t{o}\t{m}\n")
    head = [f"Estudos com fatores curados: {len(design)} | amostras excluidas: {len(seen)} "
            f"em {len({s for s, _ in seen})} estudos", ""]
    txt = "\n".join(head + L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
