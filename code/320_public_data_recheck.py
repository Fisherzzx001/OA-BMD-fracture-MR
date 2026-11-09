#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F2 / F5 / F10 public-data availability recheck (2026-11-09).

Read-only verification that the three unexecuted registered arms are still
unexecutable with *fit-for-purpose* public summary statistics.  It queries the
GWAS Catalog v2 REST API on exact ``reportedTrait`` terms, merges a previously
retrieved 442-study corpus, classifies candidate studies by registered test,
and writes an internal audit record.  No GWAS file is downloaded and no
reported beta/SE/p/OR/CI is touched.

    python3 code/320_public_data_recheck.py            # full recheck + record
    python3 code/320_public_data_recheck.py --check    # summary only

Honesty boundary: the existence of a public non-European / ancestry-stratified
GWAS is NOT equivalent to execution of the registered F10 test.  The register
requires de novo, ancestry-matched BMD instruments and a site-matched
(knee-OA) outcome; a sign-concordance probe built from 31 European instruments
against an all-site OA or heel-QUS outcome would be a portability check, not a
replication, and is deliberately not reported.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import json
import re
import socket
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INTERNAL = ROOT / "compliance_internal"
DEPOSIT = ROOT / "code" / "output" / "genetics"
# The deposited copy is the authoritative one: the GitHub/Zenodo archive
# includes code/ and code/output/ but excludes the gitignored compliance_internal/.
CORPUS = DEPOSIT / "gwas_catalog_v2_corpus_2026-11-09.json"
INTERNAL_CORPUS = INTERNAL / "gwas_catalog_v2_corpus_2026-11-09.json"
OUT_JSON = DEPOSIT / "L7_F2F5F10_public_data_recheck.json"
OUT_JSON_INTERNAL = INTERNAL / "public_data_recheck_2026-11-09.json"
OUT_MD = INTERNAL / "public_data_recheck_2026-11-09.md"
LEGACY_CORPUS = ROOT / "_tmp_gwas_hits.json"

API = "https://www.ebi.ac.uk/gwas/api/v2/studies"
HDR = {"User-Agent": "OA-BMD-fracture-MR/public-data-recheck (read-only)",
       "Accept": "application/json"}
TERMS = [
    "volumetric bone mineral density",
    "trabecular bone mineral density",
    "cortical bone mineral density",
    "bone area",
    "bone width",
    "osteophyte",
    "heel bone mineral density",
]

# ---- registered-test specifications -------------------------------------
SPECS = {
    "F2_volume": {
        "registered": "Areal BMD replaced by a volumetric / area-conditioned BMD "
                      "instrument, tested against the same hip and knee OA outcomes.",
        "need": "Common-SNP GWAS of true volumetric BMD (vBMD), cortical/trabecular "
                "vBMD, bone area or bone width with full summary statistics and "
                "enough genome-wide significant loci to build >=10 independent "
                "instruments; ideally site-matched to the OA outcomes.",
    },
    "F5_hip_osteophyte": {
        "registered": "Stratify the OA outcome by osteophyte-predominant versus "
                      "joint-space-narrowing subtype; registered endpoint is hip osteophyte.",
        "need": "Hip-osteophyte GWAS with common-SNP summary statistics (the OARSI 2024 "
                "abstract, 52,024 participants, IGFBP3 locus, has no deposited file).",
    },
    "F10_cross_ancestry": {
        "registered": "Re-run BMD->OA MR in non-European (EAS) OA with "
                      "ancestry-matched instruments / ancestry-stratified outcome.",
        "need": "De novo non-European BMD instruments plus an ancestry-stratified "
                "KNEE-OA outcome (Hatzikotoulas 2025 knee OA, GCST90566800), or "
                "an adequately powered non-European total-body DXA BMD GWAS "
                "paired with the same outcome definition.",
    },
}


def fetch(term: str, size: int = 200, timeout: int = 40, tries: int = 3) -> dict:
    """One exact reportedTrait query; returns {accession: study}.

    EBI is rate-limited from this sandbox, so a bounded retry is required for
    the record to be complete rather than silently corpus-only.
    """
    url = API + "?" + urllib.parse.urlencode({"reportedTrait": term, "size": size})
    last = None
    for attempt in range(tries):
        try:
            socket.setdefaulttimeout(timeout)
            req = urllib.request.Request(url, headers=HDR)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = json.loads(resp.read())
            studies = payload.get("_embedded", {}).get("studies", [])
            return {s["accessionId"]: s for s in studies if s.get("accessionId")}
        except Exception as exc:                            # noqa: BLE001
            last = exc
    raise RuntimeError(f"{type(last).__name__}: {last} after {tries} attempts")


def load_corpus() -> dict:
    """Prefer the deposited corpus, then the internal copy, then the session cache."""
    for p in (CORPUS, INTERNAL_CORPUS, LEGACY_CORPUS):
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))
    return {}


def study_line(acc: str, s: dict) -> dict:
    return {
        "accession": acc,
        "reported_trait": s.get("reportedTrait", ""),
        "first_author": s.get("firstAuthor", ""),
        "publication_date": s.get("publicationDate", ""),
        "ancestry": s.get("discoverySampleAncestry", []),
        "genotyping": s.get("genotypingTechnologies", []),
        "has_sumstats": str(s.get("summaryStatistics", "NA")).startswith("http"),
        "sumstats": s.get("summaryStatistics", "NA") if str(
            s.get("summaryStatistics", "NA")).startswith("http") else None,
    }


def pick(corpus: dict, predicate, limit: int = 60) -> list:
    out = [study_line(a, s) for a, s in corpus.items() if predicate(a, s)]
    return sorted(out, key=lambda x: x["accession"])[:limit]


def classify(corpus: dict) -> dict:
    def txt(s):
        return (s.get("reportedTrait", "") + " " + " ".join(
            s.get("genotypingTechnologies", []))).lower()

    def anc(s):
        return " ".join(s.get("discoverySampleAncestry") or [])

    f2 = pick(corpus, lambda a, s: re.search(
        r"volumetric|bone area|bone width|bone size|trabecular|cortical|femoral neck width|bone geometry",
        txt(s), re.I))
    f5 = pick(corpus, lambda a, s: "osteophyt" in txt(s))
    f10_bmd = pick(corpus, lambda a, s: re.search(
        r"bone mineral density|bone density|ebmd", txt(s), re.I)
        and re.search(r"asian|african|hispanic|admixed|greater middle eastern", anc(s), re.I))
    f10_oa = pick(corpus, lambda a, s: "osteoarthritis" in txt(s)
                  and re.search(r"asian|african|hispanic", anc(s), re.I))
    return {"F2_volume": f2, "F5_hip_osteophyte": f5,
            "F10_nonEUR_BMD": f10_bmd, "F10_ancestry_stratified_OA": f10_oa}


# ---- hand-checked fit-for-purpose adjudication ---------------------------
FINDINGS = [
    {
        "test": "F2_volume", "accession": "GCST90826226..GCST90826232",
        "label": "Garg 2026 AJHG CNV PheWAS: pelvis / legs / femoral-neck bone area",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Bone-area traits are measured, but the variants are copy-number "
                  "burden calls from sequencing read depth (N=39,704 EUR), not common "
                  "SNPs. They cannot serve as the areal-BMD-replacing volumetric "
                  "instrument set and would not overlap the 31 BMD instruments.",
    },
    {
        "test": "F2_volume", "accession": "GCST007014",
        "label": "Nielson 2016 trabecular lumbar-spine BMD, N=12,287 EUR",
        "public": False,
        "verdict": "NOT FIT",
        "reason": "True volumetric/trabecular BMD phenotype, but summaryStatistics=NA "
                  "in the Catalog; no instrument set can be reconstructed.",
    },
    {
        "test": "F2_volume", "accession": "GCST010299",
        "label": "Zhang 2020 'BMD (hip) and hip bone size' bivariate GWAS",
        "public": False,
        "verdict": "NOT FIT",
        "reason": "Bivariate (pleiotropy-selected) statistic with no deposited file; "
                  "even if available it is not a pure volumetric BMD instrument.",
    },
    {
        "test": "F5_hip_osteophyte", "accession": "GCST90080462/64, GCST90084448/50",
        "label": "Backman 2021 UKB exome: knee / any-site osteophyte ICD10 burden",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Exome-wide rare coding-variant / gene-burden files, not common-SNP "
                  "GWAS; 0/31 exact chr:bp matches against the locked BMD instruments "
                  "(verified 2026-11-09), and the endpoint is knee or any-site, not hip.",
    },
    {
        "test": "F5_hip_osteophyte", "accession": "GCST005278, GCST90080463, GCST90084449",
        "label": "Lumbar spondylosis osteophytes / ankle-foot osteophyte",
        "public": False,
        "verdict": "NOT FIT",
        "reason": "Wrong anatomical site and (for GCST005278) no deposited summary "
                  "statistics; the registered hip-osteophyte GWAS remains only an "
                  "OARSI 2024 conference abstract.",
    },
    {
        "test": "F10_cross_ancestry", "accession": "GCST90692328",
        "label": "Pan-UKB 2025 Central/South Asian heel BMD (UKB f.3148), N=2,633",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Non-European BMD summary statistics do exist (harmonised 484 MB, "
                  "GRCh37, .tbi), but N=2,633 heel quantitative ultrasound is a "
                  "different bone phenotype from total-body DXA areal BMD and is far "
                  "too small for de novo instrument discovery; using European "
                  "instruments against it would be a portability check, not replication.",
    },
    {
        "test": "F10_cross_ancestry", "accession": "GCST005348 / GCST007691",
        "label": "Medina-Gomez 2018 total-body BMD (incl. 1,333 AA) / GEFOS 2012 FN BMD",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Deposited files are pooled multi-ancestry or European-dominant "
                  "meta-analyses, not ancestry-specific instrument-discovery GWAS.",
    },
    {
        "test": "F10_cross_ancestry", "accession": "GCST90134280/81/82/87",
        "label": "MVP+UKB McDonald 2022 ancestry-stratified OA (SAS 7,162; EAS 1,704; AFR 7,188/35,753)",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Ancestry-stratified OA summary statistics DO exist, but the trait is "
                  "all-site osteoarthritis, not the knee-OA endpoint of the registered "
                  "test; harmonised files are 0.79-1.0 GB and no ancestry-matched BMD "
                  "instrument set exists to pair with them.",
    },
    {
        "test": "F10_cross_ancestry", "accession": "GCST90727231",
        "label": "Genes & Health 2026 British South Asian OA excluding spine, N=44,026",
        "public": True,
        "verdict": "NOT FIT",
        "reason": "Large non-European OA GWAS with deposited harmonised statistics, but "
                  "exome-wide sequencing of rare coding variants (GRCh38), an "
                  "'excluding spine' all-site endpoint, and no knee-OA definition; "
                  "common-SNP BMD instruments cannot be harmonised to it.",
    },
    {
        "test": "F10_cross_ancestry", "accession": "GCST90566811",
        "label": "Previously staged 'European-restricted knee OA' file",
        "public": True,
        "verdict": "MISIDENTIFIED",
        "reason": "This is the X-chromosome knee-OA GWAS (all rows chr23), not a "
                  "European-restricted autosomal outcome; 0/31 autosomal BMD "
                  "instruments match. The S11 sentence claiming a staged EUR file "
                  "was factually wrong and has been corrected.",
    },
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="print summary, do not write")
    args = ap.parse_args()

    corpus = load_corpus()
    fresh, errors = {}, []
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        futs = {pool.submit(fetch, t): t for t in TERMS}
        for f in cf.as_completed(futs):
            t = futs[f]
            try:
                fresh.update(f.result())
            except Exception as exc:                       # network-limited sandbox
                errors.append(f"{t}: {type(exc).__name__}: {exc}")
    merged = dict(corpus)
    merged.update({k: v for k, v in fresh.items() if k not in merged})
    new_ids = sorted(set(fresh) - set(corpus))

    classes = classify(merged)
    record = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "source": "GWAS Catalog REST v2 /studies?reportedTrait=<term> (exact); "
                  "legacy v1 /gwas/rest/api returns HTTP 410 Gone",
        "corpus_studies": len(corpus),
        "fresh_query_terms": TERMS,
        "fresh_studies": len(fresh),
        "new_accessions_since_corpus": new_ids,
        "query_errors": errors,
        "registered_specs": SPECS,
        "candidate_counts": {k: len(v) for k, v in classes.items()},
        "adjudication": FINDINGS,
        "conclusion": {
            "F2_volume": "UNEXECUTED - no fit-for-purpose common-SNP volumetric BMD "
                         "summary statistics; geometry arm and conditional MVMR stand.",
            "F5_hip_osteophyte": "UNEXECUTED - no hip-osteophyte common-SNP GWAS; "
                                 "site-stratified surrogate (Table S29) stands.",
            "F10_cross_ancestry": "UNEXECUTED - non-European / ancestry-stratified "
                                  "files now exist but none are fit for purpose; "
                                  "Table S32 ancestry-composition proxy stands and is "
                                  "NOT a replication; S23 verdict 9/10 unchanged.",
        },
        "wording_note": "Blanket phrases such as 'no public ancestry-stratified data' "
                        "and 'no non-European BMD exposure summary statistics are "
                        "available' are imprecise as of 2026-11-09 and were replaced "
                        "by 'no fit-for-purpose / no ancestry-stratified summary "
                        "statistics for this outcome study'.",
    }

    lines = [
        "# F2 / F5 / F10 public-data availability recheck (2026-11-09)",
        "",
        "Internal verification record — not part of the submission. Read-only; no "
        "reported beta/SE/p/OR/CI was modified.",
        "",
        f"Corpus: {len(corpus)} studies; fresh exact-term queries: {len(TERMS)} "
        f"({len(fresh)} studies, {len(new_ids)} new accessions). Errors: {len(errors)}.",
        f"New accessions: {', '.join(new_ids) if new_ids else 'none'}.",
        "",
        "## Registered specifications",
    ]
    for k, v in SPECS.items():
        lines += [f"- **{k}** — {v['registered']} *Need:* {v['need']}"]
    lines += ["", "## Candidate census",
              "| axis | candidate studies in corpus |", "|---|---|"]
    for k, v in classes.items():
        lines.append(f"| {k} | {len(v)} |")
    lines += ["", "## Fit-for-purpose adjudication", ""]
    for f in FINDINGS:
        lines += [
            f"### [{f['test']}] {f['accession']} — {f['verdict']}",
            f"- {f['label']}",
            f"- Deposited summary statistics: {'yes' if f['public'] else 'no / NA'}",
            f"- {f['reason']}", "",
        ]
    lines += ["## Conclusion", ""]
    for k, v in record["conclusion"].items():
        lines.append(f"- **{k}**: {v}")
    lines += ["", "## Wording", "", record["wording_note"], "",
              "## Honesty boundary", "",
              "The existence of a public non-European or ancestry-stratified GWAS is "
              "not equivalent to executing the registered F10 test. A sign-concordance "
              "probe built from 31 European instruments against an all-site OA or "
              "heel-QUS outcome would be a trans-ancestry portability check (and is "
              "underpowered: N=2,633 for the only non-European BMD file). It was "
              "considered and deliberately not run or reported; a streaming prototype "
              "(code/319) was abandoned for this reason.", ""]

    if not args.check:
        INTERNAL.mkdir(parents=True, exist_ok=True)
        DEPOSIT.mkdir(parents=True, exist_ok=True)
        if not CORPUS.exists():
            src = next((p for p in (INTERNAL_CORPUS, LEGACY_CORPUS) if p.exists()), None)
            if src:
                CORPUS.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        payload = json.dumps(record, ensure_ascii=False, indent=2)
        OUT_JSON.write_text(payload, encoding="utf-8")
        OUT_JSON_INTERNAL.write_text(payload, encoding="utf-8")
        OUT_MD.write_text("\n".join(lines), encoding="utf-8")
        print(f"[write] {OUT_JSON.relative_to(ROOT)} (deposited)")
        print(f"[write] {OUT_JSON_INTERNAL.relative_to(ROOT)} (internal)")
        print(f"[write] {OUT_MD.relative_to(ROOT)} (internal)")
    print(f"[summary] corpus={len(corpus)} fresh={len(fresh)} new={len(new_ids)} "
          f"errors={len(errors)}")
    for k, v in record["conclusion"].items():
        print(f"  {k}: {v.split(' - ')[0]}")
    if errors:
        print("[warn] queries failed (corpus still used):", "; ".join(errors[:3]),
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
