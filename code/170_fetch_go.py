#!/usr/bin/env python3
"""Fetch GO-term -> human gene-symbol memberships from QuickGO (paginated search).
Saves to /g/OA_B_PROJECT_cache/t2/go_genesets.json and .txt.
Self-contained: shells out to curl (proxy-aware) and parses JSON.
"""
import subprocess, json, os, time, sys

OUT_DIR = "/g/OA_B_PROJECT_cache/t2"
JSON_OUT = os.path.join(OUT_DIR, "go_genesets.json")
TXT_OUT = os.path.join(OUT_DIR, "go_genesets.txt")

# GO term -> human-readable label (bone vs cartilage focus)
TERMS = {
    "GO:0060348": "bone_development",
    "GO:0051216": "cartilage_development",
    "GO:0001503": "ossification",
    "GO:0002062": "chondrocyte_differentiation",
    "GO:0060347": "bone_morphogenesis",
    "GO:0030198": "extracellular_matrix_organization",
    "GO:0001501": "skeletal_system_development",
    "GO:0048705": "skeletal_joint_morphogenesis",
    "GO:0036072": "bone_remodeling",
}

def curl_json(url, timeout=60, tries=6):
    for i in range(tries):
        try:
            p = subprocess.run(
                ["curl", "-s", "--max-time", str(timeout), url],
                capture_output=True, text=True)
            if p.returncode == 0 and p.stdout.strip():
                return json.loads(p.stdout)
        except Exception as e:
            sys.stderr.write(f"  try {i+1} err {e}\n")
        time.sleep(3)
    return None

def fetch_term(term):
    genes = {}
    page = 1
    while True:
        url = (f"https://www.ebi.ac.uk/QuickGO/services/annotation/search"
               f"?goId={term}&taxonId=9606&limit=200&page={page}&format=json")
        d = curl_json(url)
        if not d:
            break
        res = d.get("results", [])
        if not res:
            break
        for r in res:
            sym = r.get("dbObjectSymbol") or r.get("symbol")
            if sym:
                genes[sym] = genes.get(sym, 0) + 1
        # page info
        pi = d.get("pageInfo", {})
        total = pi.get("total") or d.get("numberOfHits") or 0
        if len(genes) >= total and total > 0:
            break
        if len(res) < 200:
            break
        page += 1
        if page > 50:
            break
        time.sleep(0.5)
    return set(genes.keys())

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    out = {}
    for term, label in TERMS.items():
        print(f"[fetch] {term} ({label}) ...", flush=True)
        syms = fetch_term(term)
        out[term] = {"label": label, "genes": sorted(syms)}
        print(f"    -> {len(syms)} genes", flush=True)
        time.sleep(1)
    with open(JSON_OUT, "w") as f:
        json.dump(out, f, indent=1)
    with open(TXT_OUT, "w") as f:
        for term, v in out.items():
            f.write(f">{term} {v['label']} n={len(v['genes'])}\n")
            f.write("\t".join(v["genes"]) + "\n")
    print(f"SAVED {JSON_OUT} and {TXT_OUT}")

if __name__ == "__main__":
    main()
