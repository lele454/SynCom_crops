#!/usr/bin/env python3
# [S5.2] Growth-promotion traits from KOfam HMM profiles
# Supplementary 5 of the manuscript. See README.md for the full map.
"""
Triangulacao funcional, etapa 2: tracos PGP nos genomas (KOfam HMM + limiares oficiais).

1. Baixa ko_list (limiares) e extrai so os HMMs necessarios de profiles.tar.gz do KOfam
   (lido em fluxo: o arquivo de ~1,4 GB nao e gravado no disco)
2. Genomas sem proteinas: predicao de genes com pyrodigal
3. hmmsearch (pyhmmer) com o limiar de cada KO (score_type full/domain do ko_list)
4. Traco presente se a regra do traco e satisfeita (ex.: fixacao de N = nifH + nifD + nifK)
5. Por genero: fracao dos genomas com cada traco; por consorcio: tracos do centro, dos
   parceiros, uniao e tracos complementares (ausentes no centro, presentes num parceiro)

Requer: pip install --user pyhmmer pyrodigal
Saidas: <root>/genomas/pgp/ (ko_por_genoma.tsv, tracos_por_genoma.tsv, tracos_por_genero.tsv,
        tracos_por_consorcio.tsv, resumo.txt)
Uso: python3 pgp_traits.py --root ~/syncom [--cpus 4]
"""
import argparse, glob, gzip, io, os, sys, tarfile, urllib.request
import pandas as pd

KOFAM = "https://www.genome.jp/ftp/db/kofam/"
# KO -> palavra esperada na definicao do KEGG (conferencia automatica)
KOS = {
    "K02588": "nifH", "K02586": "nifD", "K02591": "nifK",
    "K01505": "1-aminocyclopropane",
    "K06136": "pqqB", "K06137": "pqqC", "K06138": "pqqD", "K06139": "pqqE",
    "K00117": "glucose dehydrogenase",
    "K01077": "phoA", "K01113": "phoD", "K01093": "appA",
    "K04103": "ipdC", "K00466": "iaaM",
    "K01183": "chitinase",
    "K10814": "hcnA", "K10815": "hcnB", "K10816": "hcnC",
    "K01575": "alsD", "K00004": "butanediol",
}
# traco: lista de grupos; cada grupo = (KOs, minimo presente)
TRAITS = {
    "fixacao_N": [(["K02588", "K02586", "K02591"], 3)],
    "ACC_desaminase": [(["K01505"], 1)],
    "PQQ_solub_fosfato": [(["K06136", "K06137", "K06138", "K06139"], 3)],
    "glicose_desidrogenase_gcd": [(["K00117"], 1)],
    "fosfatases_fitase": [(["K01077", "K01113", "K01093"], 1)],
    "AIA_via_IPyA_ou_IAM": [(["K04103", "K00466"], 1)],
    "quitinase": [(["K01183"], 1)],
    "HCN": [(["K10814", "K10815", "K10816"], 3)],
    "acetoina_butanodiol": [(["K01575"], 1), (["K00004"], 1)],
}


def read_fasta(path):
    name, seq = None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def get_hmms(d, profiles=None):
    hmm = os.path.join(d, "pgp_kos.hmm"); kl = os.path.join(d, "ko_list.tsv")
    if not os.path.exists(kl):
        print("Baixando ko_list...", flush=True)
        with urllib.request.urlopen(KOFAM + "ko_list.gz") as r:
            open(kl, "wb").write(gzip.decompress(r.read()))
    if not os.path.exists(hmm):
        print("Extraindo HMMs do KOfam em fluxo (~1,4 GB baixados, nada gravado alem dos HMMs)...", flush=True)
        want = {f"{k}.hmm" for k in KOS}; got = {}
        src = open(profiles, "rb") if profiles else urllib.request.urlopen(KOFAM + "profiles.tar.gz")
        with src as r:
            with tarfile.open(fileobj=r, mode="r|gz") as tf:
                for m in tf:
                    b = os.path.basename(m.name)
                    if b in want:
                        got[b] = tf.extractfile(m).read()
                        if len(got) == len(want):
                            break
        with open(hmm, "wb") as fo:
            for b in sorted(got):
                fo.write(got[b])
        miss = want - set(got)
        if miss:
            print("HMMs nao encontrados: " + ", ".join(sorted(miss)))
    ko = pd.read_csv(kl, sep="\t").set_index("knum")
    return hmm, ko


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--profiles", default=None, help="profiles.tar.gz do KOfam ja baixado (pula o download)")
    a = ap.parse_args()
    import pyhmmer, pyrodigal
    root = os.path.abspath(os.path.expanduser(a.root))
    gd = os.path.join(root, "genomas"); od = os.path.join(gd, "pgp"); os.makedirs(od, exist_ok=True)
    hmmf, ko = get_hmms(od, a.profiles)

    L = ["Conferencia dos KOs (definicao KEGG):"]
    for k, kw in KOS.items():
        df = str(ko.definition.get(k, "AUSENTE")) if k in ko.index else "AUSENTE"
        ok = kw.lower() in df.lower()
        L.append(f"  {k} {'ok ' if ok else 'CONFERIR'} {df[:90]}")
    print("\n".join(L), flush=True)

    # genomas
    sel = pd.read_csv(os.path.join(gd, "selecionados.tsv"), sep="\t")
    faas = {}
    for _, r in sel.iterrows():
        base = os.path.join(gd, r.genus)
        faa = glob.glob(os.path.join(base, f"{r.acc}*.faa"))
        fna = glob.glob(os.path.join(base, f"{r.acc}*.fna"))
        if faa:
            faas[(r.genus, r.acc)] = faa[0]
        elif fna:
            out = os.path.join(base, f"{r.acc}.prodigal.faa")
            if not os.path.exists(out):
                print(f"pyrodigal: {r.genus} {r.acc}", flush=True)
                contigs = [(n, s.encode()) for n, s in read_fasta(fna[0])]
                gf = pyrodigal.GeneFinder(meta=sum(len(s) for _, s in contigs) < 100000)
                if not gf.meta:
                    gf.train(*[s for _, s in contigs])
                with open(out, "w") as fo:
                    for n, s in contigs:
                        gf.find_genes(s).write_translations(fo, sequence_id=n)
            faas[(r.genus, r.acc)] = out
    print(f"Genomas com proteinas: {len(faas)}/{len(sel)}", flush=True)

    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(hmmf) as f:
        hmms = list(f)
    nm = lambda x: x.decode() if isinstance(x, bytes) else x
    rows = []
    for (g, acc), faa in faas.items():
        with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alpha) as sf:
            seqs = sf.read_block()
        found = {}
        for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=a.cpus, E=1e-5):
            q = nm(getattr(hits, "query_name", None) or hits.query.name)
            thr = ko.threshold.get(q); st = ko.score_type.get(q, "full")
            try:
                thr = float(thr)
            except (TypeError, ValueError):
                thr = None
            best = 0.0
            for h in hits:
                s = h.score if st == "full" else h.best_domain.score
                best = max(best, s)
            found[q] = best >= thr if thr is not None else best > 0
        rows.append({"genus": g, "acc": acc, **{k: int(found.get(k, False)) for k in KOS}})
    K = pd.DataFrame(rows)
    K.to_csv(os.path.join(od, "ko_por_genoma.tsv"), sep="\t", index=False)

    T = K[["genus", "acc"]].copy()
    for t, groups in TRAITS.items():
        ok = None
        for kos, mn in groups:
            v = K[kos].sum(axis=1) >= mn
            ok = v if ok is None else (ok | v)
        T[t] = ok.astype(int)
    T.to_csv(os.path.join(od, "tracos_por_genoma.tsv"), sep="\t", index=False)
    G = T.drop(columns="acc").groupby("genus").mean()
    G.insert(0, "n_genomas", T.groupby("genus").size())
    G.to_csv(os.path.join(od, "tracos_por_genero.tsv"), sep="\t", float_format="%.2f")

    sl = pd.read_csv(os.path.join(root, "merged", "meta", "consorcios", "shortlist.tsv"), sep="\t")
    tr = list(TRAITS)
    C = []
    for _, r in sl.iterrows():
        hub, parts = r.centro, r.parceiros.split(" + ")
        has = lambda g, t: g in G.index and G.loc[g, t] >= 0.5
        h = [t for t in tr if has(hub, t)]
        p = [t for t in tr if any(has(x, t) for x in parts)]
        comp = [t for t in p if t not in h]
        C.append({"id": r.id, "centro": hub, "parceiros": r.parceiros,
                  "n_tracos_centro": len(h), "n_tracos_uniao": len(set(h) | set(p)),
                  "tracos_centro": ";".join(h), "tracos_complementares": ";".join(comp),
                  "sem_genoma": ";".join(x for x in [hub] + parts if x not in G.index)})
    C = pd.DataFrame(C).sort_values(["n_tracos_uniao", "n_tracos_centro"], ascending=False)
    C.to_csv(os.path.join(od, "tracos_por_consorcio.tsv"), sep="\t", index=False)

    L += ["", "Fracao de genomas com cada traco, por genero:",
          G.to_string(float_format=lambda v: f"{v:.2f}"), "",
          "Consorcios (traco conta se >= 50% dos genomas do genero o tem):",
          C[["id", "centro", "parceiros", "n_tracos_centro", "n_tracos_uniao", "tracos_complementares"]]
          .to_string(index=False)]
    txt = "\n".join(L)
    open(os.path.join(od, "resumo.txt"), "w").write(txt)
    print("\n" + txt)


if __name__ == "__main__":
    main()
