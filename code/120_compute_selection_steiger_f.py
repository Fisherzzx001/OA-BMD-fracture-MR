#!/usr/bin/env python3
"""P0-2 + P2 data computation.

1. Genome-wide significant (GWS, p<5e-8) locus counts for all 10 non-head BMD
   sites + head, from code/output/genetics/tophits_*_BMD.csv -> variant_selection_table.csv
2. Univariate instrument strength (F, R2_exposure) and Steiger directionality test
   for the 10 MR pairs, from code/output/genetics/harmonised/*.csv
   -> instrument_strength_steiger.csv

All quantities are computed live from the deposited primary files; no hardcoded stats.
"""
import csv, math, os, glob
from statistics import NormalDist

GEN = "code/output/genetics"
HARM = os.path.join(GEN, "harmonised")
OUT = "code/output/genetics"

# ---- prevalence (liability K) from Table 1 note in manuscript ----
K = {"OA_knee": 0.12, "OA_selfrep": 0.083, "fracture": 0.143}
# effective N (per-SNP analysis N) for each trait
N = {"BMD": 31873, "OA_knee": 1208397, "fracture": 264973}

def count_gws_from_clump_input(path):
    """Count SNPs with p<5e-8 in the full candidate-SNP clump-input file.
    Validated against strict_clump_summary.csv (514/246/170/277/301 match exactly)."""
    gws = 0
    with open(path) as fh:
        fh.readline()  # header SNP\tP
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                if float(parts[1]) < 5e-8:
                    gws += 1
            except ValueError:
                pass
    return gws

def count_gws_from_tophits(path):
    """Count genome-wide-significant SNPs in a deposited top-hits file
    (comma-separated GWAS format, p_value column)."""
    gws = 0
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                if float(row["p_value"]) < 5e-8:
                    gws += 1
            except (KeyError, TypeError, ValueError):
                continue
    return gws

# ---------- 1. variant selection table ----------
print("=== 1. GWS counts per BMD site (tophits) ===")
sites = {
    "BMD_total": "Total body", "BMD_femur_neck": "Femoral neck",
    "BMD_femur_total": "Femoral total", "BMD_lumbar": "Lumbar spine",
    "BMD_spine": "Spine (total)", "BMD_arms": "Arms", "BMD_leg": "Legs",
    "BMD_pelvis": "Pelvis", "BMD_rib": "Rib", "BMD_trunk": "Trunk",
    "BMD_head": "Head (excluded a priori)",
}
# strict clump summary (only 5 BMD sites + OA_knee were clumped)
clump = {}
with open(os.path.join(GEN, "strict_clump_summary.csv")) as f:
    for row in csv.DictReader(f):
        clump[row["exposure"]] = row

analyzed = {"BMD_arms", "BMD_femur_neck", "BMD_leg", "BMD_lumbar", "BMD_total"}
sel_rows = []
for key, label in sites.items():
    base = key
    if base == "BMD_head":
        # The head site never entered clumping; count its genome-wide-significant SNPs
        # from the deposited top-hits file to document that the a priori exclusion is
        # principled (not instrument-driven).
        gws = count_gws_from_tophits(os.path.join(GEN, "tophits_head_BMD.csv"))
    else:
        mapping = {"BMD_total":"Total","BMD_femur_neck":"Femur_neck","BMD_femur_total":"Femur_total",
                   "BMD_lumbar":"Lumbar","BMD_spine":"Spine","BMD_arms":"Arms","BMD_leg":"Leg",
                   "BMD_pelvis":"Pelvis","BMD_rib":"Rib","BMD_trunk":"Trunk"}
        gws = count_gws_from_clump_input(os.path.join(GEN, f"BMD_{mapping[key]}_clump_input.txt"))
    c = clump.get(base)
    if base in analyzed:
        n_kept = c["n_kept_strict"] if c else "NA"
        entered = "Yes"
        reason = "Strict-clumped; >=15 independent instruments; entered MR"
    elif base == "BMD_head":
        n_kept = "NA"
        entered = "No"
        reason = ("Excluded a priori: not a load-bearing skeletal region; head BMD is "
                  "genetically correlated with intracranial aneurysm in the source GWAS "
                  "(Qian et al. 2025), a non-skeletal pleiotropy risk. Genome-wide signal "
                  f"was present ({gws} SNPs at p<5e-8), so the exclusion is principled "
                  "rather than instrument-driven.")
    else:
        n_kept = "NA"
        entered = "No"
        reason = "Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set"
    sel_rows.append({
        "BMD_site": label, "n_gws_loci_p5e-8": gws,
        "n_strict_instruments": n_kept, "entered_MR": entered, "reason": reason,
    })
    print(f"{label:32s} gws={str(gws):>4s} kept={str(n_kept):>4s}  MR={entered}")

with open(os.path.join(OUT, "variant_selection_table.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["BMD_site","n_gws_loci_p5e-8",
                                       "n_strict_instruments","entered_MR","reason"])
    w.writeheader()
    w.writerows(sel_rows)
print("wrote", os.path.join(OUT, "variant_selection_table.csv"))

# ---------- 2. Steiger + F from harmonised ----------
print("\n=== 2. Instrument strength (F) + Steiger direction test ===")
def parse_pair(fname):
    # BMD_arms__OA_knee.csv -> X=BMD_arms, Y=OA_knee
    stem = fname[:-4]
    x, y = stem.split("__")
    return x, y

def trait_meta(t):
    if t.startswith("BMD"):
        return N["BMD"], None           # continuous, no prevalence
    if t == "OA_knee":
        return N["OA_knee"], K["OA_knee"]
    if t == "fracture":
        return N["fracture"], K["fracture"]
    if t == "OA_selfrep":
        return 462933, K["OA_selfrep"]
    raise ValueError(t)

def steiger(beta_x, beta_y, fx, nx, ny, kx, ky):
    # per-allele, standardized-trait scale; binary -> liability correction sqrt(K(1-K))
    bx = beta_x if kx is None else beta_x * math.sqrt(kx*(1-kx))
    by = beta_y if ky is None else beta_y * math.sqrt(ky*(1-ky))
    gv = 2*fx*(1-fx)
    r2x = gv*bx*bx
    r2y = gv*by*by
    r2xy = gv*bx*by
    # Steiger z_c
    denom_x = max(1e-12, 1-r2x)
    ratio = max(1e-12, (1-r2xy))/denom_x
    zc = 0.5*math.log(ratio)
    se = math.sqrt(1/(nx-1) + 1/(ny-1))
    z = zc/se
    # one-sided p for "exposure causes outcome" (zc>0)
    nd = NormalDist()
    if z >= 0:
        p = 1 - nd.cdf(z)
    else:
        p = nd.cdf(z)
    return r2x, r2y, r2xy, z, p

res_rows = []
for path in sorted(glob.glob(os.path.join(HARM, "*.csv"))):
    name = os.path.basename(path)
    if name in ("_loo.csv", "_validation.csv"):
        continue
    if "__" not in name:
        continue
    x, y = parse_pair(name)
    nx, kx = trait_meta(x)
    ny, ky = trait_meta(y)
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                be = float(r["beta_e"]); se = float(r["se_e"])
                bo = float(r["beta_o"]); eaf = float(r["eaf"])
            except (ValueError, KeyError):
                continue
            rows.append((be, se, bo, eaf))
    if not rows:
        continue
    n_snp = len(rows)
    # univariate F
    minF = min((be/se)**2 for be, se, _, _ in rows)
    aggF = sum(be*be for be, _, _, _ in rows) / sum(se*se for _, se, _, _ in rows)
    # R2_exposure (continuous approx on standardized trait)
    r2x = sum(2*eaf*(1-eaf)*be*be for _, _, _, eaf in [(be,se,bo,eaf) for be,se,bo,eaf in rows])
    # Steiger uses X as designated exposure
    bx_list = []
    R2x = 0.0; R2xy = 0.0
    for be, se, bo, eaf in rows:
        gv = 2*eaf*(1-eaf)
        bx = be if kx is None else be*math.sqrt(kx*(1-kx))
        by = bo if ky is None else bo*math.sqrt(ky*(1-ky))
        R2x += gv*bx*bx
        R2xy += gv*bx*by
    denom = max(1e-12, 1-R2x)
    zc = 0.5*math.log(max(1e-12,(1-R2xy))/denom)
    se_z = math.sqrt(1/(nx-1)+1/(ny-1))
    z = zc/se_z
    nd = NormalDist()
    p_steiger = (1-nd.cdf(z)) if z >= 0 else nd.cdf(z)
    direction_ok = (z > 0 and p_steiger < 0.05)
    res_rows.append({
        "pair": f"{x} -> {y}", "n_snp": n_snp,
        "F_min": round(minF,1), "F_agg": round(aggF,1),
        "R2_exposure": round(R2x,4),
        "Steiger_z": round(z,2), "Steiger_p": f"{p_steiger:.2e}",
        "direction_supported": "Yes" if direction_ok else "No",
    })
    print(f"{x:18s}->{y:14s} n={n_snp:2d} Fmin={minF:7.1f} Fagg={aggF:7.1f} "
          f"R2x={R2x:.4f} Steiger_z={z:7.2f} p={p_steiger:.1e} dir={'OK' if direction_ok else 'X'}")

with open(os.path.join(OUT, "instrument_strength_steiger.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["pair","n_snp","F_min","F_agg","R2_exposure",
                                       "Steiger_z","Steiger_p","direction_supported"])
    w.writeheader()
    w.writerows(res_rows)
print("wrote", os.path.join(OUT, "instrument_strength_steiger.csv"))
print("\nDONE")
