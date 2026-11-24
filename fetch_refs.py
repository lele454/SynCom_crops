#!/usr/bin/env python3
# [S1.16] Publication linked to each deposit
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Referencia do trabalho original de cada estudo reaproveitado.

Duas fontes, nesta ordem:
  1. ENA browser API: o vinculo PUBMED registrado no proprio deposito (dado curado)
  2. Europe PMC: busca pela acessao no texto dos artigos (casamento de texto, pode errar)

A coluna ref_source diz de onde veio cada referencia, de modo que as obtidas por
casamento de texto possam ser conferidas uma a uma. Estudos sem publicacao vinculada
ficam em branco: muitos depositos realmente nao tem artigo associado, e inventar
referencia e pior que deixar vazio.

Retomavel: roda de novo e so consulta os estudos que ainda nao estao na saida.

O Europe PMC pede um contato no User-Agent para uso automatizado. Para informa-lo:
  export SYNCOM_CONTACT="seu-email"
Sem a variavel, o script roda assim mesmo, com User-Agent generico.

Uso: python3 fetch_refs.py --tabela Tabelas/TableS1_studies.tsv --out Tabelas/study_refs.tsv
"""
import argparse, csv, difflib, json, os, re, sys, time
import urllib.error, urllib.parse, urllib.request

ENA = "https://www.ebi.ac.uk/ena/browser/api/xml/{}"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}.fcgi"
_apikey = os.environ.get("NCBI_API_KEY", "").strip()
_contato = os.environ.get("SYNCOM_CONTACT", "").strip()
UA = "syncom-meta-analysis/1.0" + (f" ({_contato})" if _contato else "")
CAMPOS = ["study", "pmid", "doi", "ref_title", "ref_journal", "ref_year", "ref_source"]


def pega(url, params=None, tentativas=3):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for i in range(tentativas):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(1.5 * (i + 1))
        except Exception:
            time.sleep(1.5 * (i + 1))
    return None


def pmid_do_ena(acc):
    """PMID registrado como referencia cruzada do estudo no ENA."""
    xml = pega(ENA.format(acc))
    if not xml:
        return None
    # <XREF_LINK><DB>PUBMED</DB><ID>12345678</ID></XREF_LINK>
    for m in re.finditer(r"<XREF_LINK>\s*<DB>([^<]+)</DB>\s*<ID>([^<]+)</ID>", xml, re.I):
        if m.group(1).strip().upper() in ("PUBMED", "PUBMED CENTRAL"):
            ids = [x.strip() for x in m.group(2).split(",") if x.strip().isdigit()]
            if ids:
                return ids[0]
    return None


def pmid_do_ncbi(acc):
    """PMID vinculado ao BioProject no NCBI. E a fonte dos depositos PRJNA, que o
    espelho do ENA nao traz."""
    p = {"db": "bioproject", "term": f"{acc}[Project Accession]", "retmode": "json"}
    if _apikey:
        p["api_key"] = _apikey
    txt = pega(EUTILS.format("esearch"), p)
    if not txt:
        return None
    try:
        ids = json.loads(txt).get("esearchresult", {}).get("idlist", []) or []
    except Exception:
        return None
    if not ids:
        return None
    p = {"dbfrom": "bioproject", "db": "pubmed", "id": ids[0], "retmode": "json"}
    if _apikey:
        p["api_key"] = _apikey
    txt = pega(EUTILS.format("elink"), p)
    if not txt:
        return None
    try:
        for ls in json.loads(txt).get("linksets", []):
            for db in ls.get("linksetdbs", []) or []:
                for i in db.get("links", []) or []:
                    if str(i).isdigit():
                        return str(i)
    except Exception:
        return None
    return None


def epmc(query, n=5):
    txt = pega(EPMC, {"query": query, "format": "json", "pageSize": n, "resultType": "core"})
    if not txt:
        return []
    try:
        return json.loads(txt).get("resultList", {}).get("result", []) or []
    except Exception:
        return []


def registro(r, fonte):
    return {"pmid": r.get("pmid", ""), "doi": r.get("doi", ""),
            "ref_title": (r.get("title") or "").rstrip("."),
            "ref_journal": ((r.get("journalInfo") or {}).get("journal") or {}).get("title", ""),
            "ref_year": r.get("pubYear", ""), "ref_source": fonte}


def normaliza(s):
    return re.sub(r"[^a-z0-9 ]+", " ", str(s).lower()).split()


def por_titulo(titulo, minimo):
    """Casa o titulo do deposito com o titulo do artigo. Alcanca tambem artigos pagos,
    que o Europe PMC indexa por titulo mas nao por texto completo."""
    pal = normaliza(titulo)
    if len(pal) < 4:
        return None
    consulta = " ".join(pal[:14])
    melhor, escore = None, 0.0
    for r in epmc(f'TITLE:"{consulta}"', n=10) or epmc(" ".join(f"TITLE:{p}" for p in pal[:10]), n=10):
        t = " ".join(normaliza(r.get("title") or ""))
        s = difflib.SequenceMatcher(None, " ".join(pal), t).ratio()
        if s > escore:
            melhor, escore = r, s
    if melhor is not None and escore >= minimo:
        return registro(melhor, f"Europe PMC title match (score {escore:.2f})")
    return None


def busca(acc, titulo="", minimo=0.60):
    for fonte, fn in (("NCBI BioProject link", pmid_do_ncbi), ("ENA cross-reference", pmid_do_ena)):
        pmid = fn(acc)
        if not pmid:
            continue
        hits = epmc(f"EXT_ID:{pmid} AND SRC:MED", n=1)
        if hits:
            return registro(hits[0], fonte)
        return {"pmid": pmid, "doi": "", "ref_title": "", "ref_journal": "", "ref_year": "",
                "ref_source": fonte}
    hits = epmc(f'"{acc}"', n=5)
    if len(hits) == 1:
        return registro(hits[0], "Europe PMC text match")
    if len(hits) > 1:
        return registro(hits[0], f"Europe PMC text match ({len(hits)} hits, check)")
    r = por_titulo(titulo, minimo)
    if r is not None:
        return r
    return {"pmid": "", "doi": "", "ref_title": "", "ref_journal": "", "ref_year": "",
            "ref_source": "no linked publication found"}


def grava(caminho, alvo, feitos):
    with open(caminho, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CAMPOS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for s in alvo:
            if s in feitos:
                w.writerow(feitos[s])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tabela", default=os.path.join("Tabelas", "TableS1_studies.tsv"))
    ap.add_argument("--out", default=os.path.join("Tabelas", "study_refs.tsv"))
    ap.add_argument("--pausa", type=float, default=0.34)
    ap.add_argument("--minimo", type=float, default=0.60,
                    help="semelhanca minima entre o titulo do deposito e o do artigo")
    ap.add_argument("--refazer-fracos", action="store_true",
                    help="reconsulta os estudos sem publicacao ou resolvidos por casamento de texto")
    a = ap.parse_args()

    with open(a.tabela, encoding="utf-8") as fh:
        linhas = list(csv.DictReader(fh, delimiter="\t"))
    usados = [r for r in linhas if r.get("status", "") != "screened out"]
    alvo = sorted(dict.fromkeys(r["study"] for r in usados))
    titulos = {r["study"]: r.get("study_title", "") for r in usados}

    feitos = {}
    if os.path.exists(a.out):
        with open(a.out, encoding="utf-8") as fh:
            feitos = {r["study"]: r for r in csv.DictReader(fh, delimiter="\t")}
    if a.refazer_fracos:
        fracos = [s for s, r in feitos.items()
                  if not r.get("ref_source", "").startswith(("NCBI", "ENA"))]
        for s in fracos:
            del feitos[s]
        print(f"reconsultando {len(fracos)} estudos sem vinculo curado")
    falta = [s for s in alvo if s not in feitos]
    print(f"{len(alvo)} estudos usados | {len(feitos)} ja consultados | {len(falta)} a consultar")

    for i, acc in enumerate(falta, 1):
        r = busca(acc, titulos.get(acc, ""), a.minimo)
        r["study"] = acc
        feitos[acc] = r
        marca = {"N": "NCBI", "E": "ENA", "n": "-"}.get(r["ref_source"][:1], "PMC")
        print(f"  [{i}/{len(falta)}] {acc:16s} {marca:4s} {r.get('doi') or r.get('pmid') or ''}")
        sys.stdout.flush()
        grava(a.out, alvo, feitos)
        time.sleep(a.pausa)

    cur = sum(1 for r in feitos.values() if r["ref_source"].startswith(("NCBI", "ENA")))
    n_pmc = sum(1 for r in feitos.values() if r["ref_source"].startswith("Europe"))
    n_chk = sum(1 for r in feitos.values() if "check" in r["ref_source"])
    n_doi = sum(1 for r in feitos.values() if r.get("doi"))
    print(f"\nvinculo curado (NCBI/ENA): {cur} | casamento de texto: {n_pmc} (conferir: {n_chk}) "
          f"| sem publicacao: {len(feitos) - cur - n_pmc}")
    print(f"com DOI: {n_doi} de {len(feitos)}")
    print(f"saida: {a.out}")


if __name__ == "__main__":
    main()
