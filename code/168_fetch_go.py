#!/usr/bin/env python3
# 168_fetch_go.py — fetch GO-term -> human gene-symbol memberships via QuickGO REST
# Used for F3 competitive gene-set test (bone-development vs cartilage-development).
import json, time, urllib.request, urllib.parse, os, sys

OUT = "/g/OA_B_PROJECT_cache/t2/go_members"
os.makedirs(OUT, exist_ok=True)

# GO terms: bone/cartilage developmental contrast + sensitivities
TERMS = {
    "bone_dev":            "GO:0060348",  # bone development
    "cartilage_dev":       "GO:0051216",  # cartilage development
    "skeletal_sys_dev":    "GO:0001501",  # skeletal system development
    "osteoblast_diff":     "GO:0045667",  # osteoblast differentiation
    "chondrocyte_diff":    "GO:0002062",  # chondrocyte differentiation
    "endochondral_oss":    "GO:0001958",  # endochondral ossification
    "cartilage_condens":   "GO:0001502",  # cartilage condensation
    "osteoblast_dev":      "GO:0001649",  # osteoblast development
}

def fetch(term, include_descendants=True):
    base = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    genes = {}
    page = 1
    while True:
        q = urllib.parse.urlencode({
            "goId": term,
            "taxonId": 9606,
            "limit": 200,
            "page": page,
            "format": "json",
            "includeDescendants": str(include_descendants).lower(),
        })
        url = base + "?" + q
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                data = json.loads(r.read().decode())
        except Exception as e:
            print("  fetch error %s: %s" % (term, e)); break
        results = data.get("results", [])
        for an in results:
            sym = an.get("symbol") or an.get("geneName") or an.get("dbObjectSymbol")
            if sym:
                genes[sym.upper()] = genes.get(sym.upper(), 0) + 1
        total = data.get("total", len(results))
        print("    %s page=%d got=%d total=%d uniq=%d" % (term, page, len(results), total, len(genes)))
        if len(results) < 200 or not results:
            break
        page += 1
        if page > 25:
            break
        time.sleep(0.4)
    return genes

if __name__ == "__main__":
    for name, go in TERMS.items():
        print("Fetching", name, go)
        g = fetch(go)
        with open(os.path.join(OUT, "%s.txt" % name), "w") as f:
            for s in sorted(g):
                f.write(s + "\n")
        print("  -> %d genes saved to %s.txt" % (len(g), name))
        time.sleep(0.6)
    print("DONE")
