#!/usr/bin/env python3
# [S1.13] Recovery of design factors from metadata
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Curadoria automatica dos fatores de desenho dos estudos que ficaram sem fatores (k = 0).

Para cada estudo sem fatores em merged/effects/study_design.tsv:
  A. atributos do BioSample com 2-20 niveis e cobertura >= 80% que o filtro automatico descartou
     (chave sem padrao reconhecido, ou bloqueada pela lista SKIP, como title/sample_name)
  B. nomes das amostras (sample_name, title, alias, library_name...): divididos em tokens por
     _ - . espaco; cada posicao com 2-20 niveis vira candidato a fator; posicoes do tipo replica
     (R1, rep2, 1, 2, 3...) sao marcadas e nao entram
  C. dossie por estudo (atributos, tokens, descricao do BioProject e resumo do artigo,
     se papers_evidence.jsonl existir) para conferencia

Saidas (<root>/merged/curadoria_auto/):
  design_auto.tsv   study, fatores (para effect_sizes.py --design)
  attrs_auto/       copia dos atributos + fatores derivados dos nomes (para effect_sizes.py --attrs)
  dossie.txt        resumo legivel por estudo
Uso: LD_PRELOAD=... python3 auto_design.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs \
        [--papers papers_evidence.jsonl]
"""
import argparse, csv, gzip, json, os, re, shutil, sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import load_attrs, detect_design, DESIGN_GROUPS, SKIP

NAME_KEYS = re.compile(r"sample.?name|^title$|alias|library.?name|sample.?title|^name$|description", re.I)
REPL = re.compile(r"^(r|rep|replicate|bio|b|s|sample|n)?[-_]?\d{1,2}[a-z]?$", re.I)
NOISE = re.compile(r"accession|ena-|insdc|checklist|submission|^_|broker|last.?update|first.?public|"
                   r"md5|ftp|bytes|read.?count|base.?count|spots|instrument|platform|library.?(strategy|source|"
                   r"selection|layout)|primer|target.?gene|seq.?meth|env.?broad|env.?local|lat.?lon$", re.I)


def samples_of(md, s, run2bs):
    with gzip.open(os.path.join(md, "genus", f"{s}.tsv.gz"), "rt") as fh:
        fh.readline()
        runs = [l.split("\t", 1)[0] for l in fh]
    return runs, [run2bs.get(r, r) for r in runs]


def tokens(v):
    return [t for t in re.split(r"[_\-\.\s]+", v.strip()) if t]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--papers", default=None)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); out = os.path.join(md, "curadoria_auto")
    ad = os.path.join(out, "attrs_auto"); os.makedirs(ad, exist_ok=True)

    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    des = list(csv.DictReader(open(os.path.join(md, "effects", "study_design.tsv")), delimiter="\t"))
    todo = [r["study"] for r in des if int(r["k"]) == 0]
    papers = {}
    if a.papers and os.path.exists(a.papers):
        for l in open(a.papers):
            if l.strip():
                e = json.loads(l); papers[e["study_accession"]] = e

    # copia todos os atributos (estudos com fatores continuam iguais)
    for f in os.listdir(a.attrs):
        shutil.copy(os.path.join(a.attrs, f), os.path.join(ad, f))

    rows_design, L = [], []
    n_ok = 0
    for s in todo:
        runs, bss = samples_of(md, s, run2bs)
        n = len(runs)
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        keys = Counter(k for b in bss for k in attrs.get(b, {}))
        cand_attr, info = [], []
        for k, c in keys.most_common():
            if NOISE.search(k):
                continue
            vals = [attrs.get(b, {}).get(k, "NA") for b in bss]
            lv = Counter(vals)
            cov = c / n
            motivo = []
            if cov < 0.8: motivo.append(f"cobertura {cov:.0%}")
            if not (2 <= len(lv) <= 20): motivo.append(f"{len(lv)} niveis")
            if SKIP.search(k): motivo.append("lista SKIP")
            elif not any(re.search(p, k, re.I) for _, p in DESIGN_GROUPS): motivo.append("sem padrao")
            ok = cov >= 0.8 and 2 <= len(lv) <= 20 and not NAME_KEYS.search(k) and not SKIP.search(k)
            if ok:
                cand_attr.append(k)
            if 2 <= len(lv) <= max(20, n):
                top = ", ".join(f"{v}({x})" for v, x in lv.most_common(6))
                info.append(f"    {'*' if ok else ' '} {k}: {len(lv)} niveis [{'; '.join(motivo) or 'ok'}] {top}")
        # tokens dos nomes
        new_rows, cand_tok, tok_info = [], [], []
        for k in [k for k in keys if NAME_KEYS.search(k)]:
            vals = [attrs.get(b, {}).get(k, "") for b in bss]
            tk = [tokens(v) for v in vals]
            if not tk:
                continue
            L_ = Counter(len(t) for t in tk).most_common(1)[0][0]
            if L_ < 2 or sum(len(t) == L_ for t in tk) < 0.8 * n:
                continue
            for p in range(L_):
                col = [t[p] if len(t) == L_ else "NA" for t in tk]
                lv = Counter(col)
                if not (2 <= len(lv) <= 20):
                    continue
                rep = all(REPL.match(x) for x in lv if x != "NA")
                name = f"auto_nome{p + 1}"
                tok_info.append(f"    {' ' if rep else '*'} {k} posicao {p + 1}: {len(lv)} niveis"
                                f"{' (replica, ignorado)' if rep else ''} {', '.join(f'{v}({x})' for v, x in lv.most_common(6))}")
                if rep:
                    continue
                cand_tok.append(name)
                for b, v in zip(bss, col):
                    new_rows.append((b, name, v))
            if cand_tok:
                break                                   # uma chave de nome basta
        if new_rows:
            with open(os.path.join(ad, f"{s}.tsv"), "a") as fo:
                if os.path.getsize(os.path.join(ad, f"{s}.tsv")) == 0:
                    fo.write("sample\tkey\tvalue\n")
                for b, k, v in sorted(set(new_rows)):
                    fo.write(f"{b}\t{k}\t{v}\n")
        fat = cand_attr + cand_tok
        if fat:
            n_ok += 1
            rows_design.append({"study": s, "fatores": ";".join(fat)})
        # dossie
        e = papers.get(s, {})
        L.append(f"=== {s} | n={n} | fatores propostos: {', '.join(fat) or 'NENHUM'}")
        L.append("  Atributos (* = proposto):"); L += info[:15] or ["    (nenhum com 2+ niveis)"]
        if tok_info:
            L.append("  Tokens dos nomes (* = proposto):"); L += tok_info
        if e.get("bioproject_desc"):
            L.append("  BioProject: " + e["bioproject_desc"][:500].replace("\n", " "))
        for p in e.get("artigos", [])[:1]:
            L.append(f"  Artigo: {p.get('title', '')} ({p.get('year', '')}) doi:{p.get('doi', '')}")
            L.append("  Resumo: " + p.get("abstract", "")[:900].replace("\n", " "))
            if p.get("metodos"):
                L.append("  Metodos: " + p["metodos"][:700].replace("\n", " "))
        L.append("")
    with open(os.path.join(out, "design_auto.tsv"), "w") as fo:
        fo.write("study\tfatores\n")
        for r in rows_design:
            fo.write(f"{r['study']}\t{r['fatores']}\n")
    head = [f"Estudos sem fatores: {len(todo)} | com fatores propostos: {n_ok} | sem proposta: {len(todo) - n_ok}",
            f"Artigos disponiveis no dossie: {sum(1 for s in todo if papers.get(s, {}).get('artigos'))}", ""]
    open(os.path.join(out, "dossie.txt"), "w").write("\n".join(head + L))
    print("\n".join(head))
    print(f"Dossie: {os.path.join(out, 'dossie.txt')}")


if __name__ == "__main__":
    main()
