#!/usr/bin/env python
# L7 F1/F2 MR analysis: BMD (exposure) -> hip-OA (F1) and knee-OA (validation/independent replication)
# Orientation: recover BMD exposure effect allele per SNP via knee-OA (== harmonised beta_o reference),
# then apply identical orientation to hip-OA.
import csv, math, json

HARM = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics/harmonised/BMD_total__OA_knee.csv"
HIP = "G:/OA_B_PROJECT_cache/l7/hip_matched.tsv"
KNEE = "G:/OA_B_PROJECT_cache/l7/knee_matched.tsv"
OUT = "G:/OA_B_PROJECT_cache/l7/l7_results.json"

# ---- load harmonised instruments (exposure BMD + reference outcome OA_knee) ----
inst = {}
with open(HARM) as f:
    r = csv.DictReader(f)
    for row in r:
        inst[row["SNP"]] = dict(
            chr=row["CHR"], bp=row["BP"],
            beta_e=float(row["beta_e"]), se_e=float(row["se_e"]), p_e=float(row["p_e"]),
            beta_o=float(row["beta_o"]), se_o=float(row["se_o"]), p_o=float(row["p_o"]),
            eaf=float(row["eaf"]), flip=float(row["flip"]),
        )

def parse_matched(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            rs = p[0]
            # columns: rs, chr, bp, E, O, OR, SE(#NA), eaf, p, variant_id, ci_u, ci_l, NCASE, NCTRL
            E, O = p[3], p[4]
            OR = float(p[5]); cu = float(p[10]); cl = float(p[11])
            ncase = float(p[12]); nctrl = float(p[13])
            # SE of ln(OR) from CI (95%): (ln(cu)-ln(cl))/3.92
            se_ln = (math.log(cu) - math.log(cl)) / 3.92
            ln_or = math.log(OR)
            d[rs] = dict(E=E, O=O, OR=OR, ln_or=ln_or, se_ln=se_ln,
                         eaf=float(p[7]), p=float(p[8]), ncase=ncase, nctrl=nctrl,
                         N=ncase+nctrl)
    return d

hip = parse_matched(HIP)
knee = parse_matched(KNEE)

# ---- Step 1: recover BMD exposure effect allele per SNP from knee ----
# harmonised beta_o is on the BMD exposure effect allele.
# knee ln(OR) is on knee E allele. If sign matches beta_o -> exposure allele = knee E; else = knee O.
orient = {}  # rs -> 'E' or 'O' meaning which knee allele is the exposure effect allele
recon_knee = {}  # rs -> reconstructed ln(OR) on exposure allele
validation = []
for rs, h in inst.items():
    if rs not in knee:
        continue
    k = knee[rs]
    sign_knee = math.copysign(1.0, k["ln_or"])
    sign_ref = math.copysign(1.0, h["beta_o"])
    if sign_knee == sign_ref:
        orient[rs] = k["E"]      # exposure effect allele == knee E
        recon = k["ln_or"]
    else:
        orient[rs] = k["O"]      # exposure effect allele == knee O (flip)
        recon = -k["ln_or"]
    recon_knee[rs] = recon
    validation.append((rs, h["beta_o"], recon, h["se_o"], k["se_ln"]))

# validation: agreement between harmonised beta_o and reconstructed knee ln(OR)
import statistics
diffs = [(v[1]-v[2]) for v in validation]
print("=== KNEE VALIDATION (reconstruct beta_o from GCST90566800) ===")
print(f"n={len(validation)} mean diff={statistics.mean(diffs):.5f} max|diff|={max(abs(d) for d in diffs):.5f}")
# correlation
xs=[v[1] for v in validation]; ys=[v[2] for v in validation]
mx=statistics.mean(xs); my=statistics.mean(ys)
cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
sx=math.sqrt(sum((x-mx)**2 for x in xs)); sy=math.sqrt(sum((y-my)**2 for y in ys))
r=cov/(sx*sy)
print(f"Pearson r between harmonised beta_o and reconstructed = {r:.4f}")

# ---- Step 2: orient hip-OA to exposure allele using recovered orient ----
hip_oriented = {}
for rs, h in inst.items():
    if rs not in hip:
        continue
    hp = hip[rs]
    exp_allele = orient.get(rs)
    if exp_allele is None:
        continue
    if hp["E"] == exp_allele:
        ln_or = hp["ln_or"]; se = hp["se_ln"]
    else:
        ln_or = -hp["ln_or"]; se = hp["se_ln"]
    hip_oriented[rs] = dict(beta_o=ln_or, se_o=se, beta_e=h["beta_e"], se_e=h["se_e"],
                            N=hp["N"], OR=math.exp(ln_or), eaf=h["eaf"])

# ---- Step 3: MR estimation (IVW fixed, IVW random, MR-Egger, weighted median) ----
def ivw_fixed(dat):
    num=0.0; den=0.0
    for d in dat:
        w = 1.0/(d["se_o"]**2)
        num += d["beta_o"]*d["beta_e"]*w
        den += d["beta_e"]**2 * w
    b=num/den; se=1.0/math.sqrt(den)
    return b, se

def ivw_random(dat):
    # DerSimonian-Laird tau^2
    num=0.0; den=0.0
    for d in dat:
        w=1.0/(d["se_o"]**2); num+=d["beta_o"]*d["beta_e"]*w; den+=d["beta_e"]**2*w
    b=num/den
    q=0.0
    for d in dat:
        w=1.0/(d["se_o"]**2)
        q += w*(d["beta_o"] - b*d["beta_e"])**2
    k=len(dat)
    if k>1:
        c=k-1  # simplified (no weight correction)
        tau2=max(0.0,(q-c)/((den) - sum((d["beta_e"]* (1.0/(d["se_o"]**2)))**2 for d in dat)/den)) if den>0 else 0.0
    else:
        tau2=0.0
    den2=sum(d["beta_e"]**2 * 1.0/(d["se_o"]**2+tau2) for d in dat)
    num2=sum(d["beta_o"]*d["beta_e"]*1.0/(d["se_o"]**2+tau2) for d in dat)
    b2=num2/den2; se2=1.0/math.sqrt(den2)
    return b2, se2, tau2

def mr_egger(dat):
    # MR-Egger: regress beta_o on beta_e with intercept
    xs=[d["beta_e"] for d in dat]; ys=[d["beta_o"] for d in dat]; ws=[1.0/(d["se_o"]**2) for d in dat]
    sw=sum(ws); swx=sum(w*x for w,x in zip(ws,xs)); swxx=sum(w*x*x for w,x in zip(ws,xs))
    swy=sum(w*y for w,y in zip(ws,ys)); swxy=sum(w*x*y for w,x,y in zip(ws,xs,ys))
    denom=sw*swxx-swx*swx
    beta_eg=(sw*swxy-swx*swy)/denom
    int_eg=(swxx*swy-swx*swxy)/denom
    se_beta=math.sqrt(sw/denom)
    se_int=math.sqrt(swxx/denom)
    return beta_eg, se_beta, int_eg, se_int

def weighted_median(dat):
    items=sorted(dat, key=lambda d: (d["beta_o"]/d["beta_e"]))
    wsum=sum(1.0/(d["se_o"]**2) for d in items)
    cum=0.0
    for d in items:
        cum+=1.0/(d["se_o"]**2)
        if cum>=wsum/2:
            return d["beta_o"]/d["beta_e"]
    return items[-1]["beta_o"]/items[-1]["beta_e"]

dat = list(hip_oriented.values())
b_ivw, se_ivw = ivw_fixed(dat)
b_ivw_r, se_ivw_r, tau2 = ivw_random(dat)
b_eg, se_eg, int_eg, se_int = mr_egger(dat)
b_wm = weighted_median(dat)

print("\n=== F1: BMD (genetic) -> hip-OA (GCST90566798) ===")
print(f"n instruments = {len(dat)}")
print(f"IVW (fixed)   : beta={b_ivw:.5f} se={se_ivw:.5f} OR/SD={math.exp(b_ivw):.4f} [{math.exp(b_ivw-1.96*se_ivw):.4f},{math.exp(b_ivw+1.96*se_ivw):.4f}] p={2*(1-0.5*(1+math.erf(b_ivw/(se_ivw*math.sqrt(2))))):.4g}")
print(f"IVW (random)  : beta={b_ivw_r:.5f} se={se_ivw_r:.5f} OR/SD={math.exp(b_ivw_r):.4f} tau2={tau2:.5f}")
print(f"MR-Egger slope: beta={b_eg:.5f} se={se_eg:.5f} OR/SD={math.exp(b_eg):.4f} intercept={int_eg:.5f} (se={se_int:.5f})")
print(f"WeightedMedian: beta={b_wm:.5f} OR/SD={math.exp(b_wm):.4f}")

# ---- also compute knee-OA MR as independent replication + validation ----
dat_k = []
for rs, h in inst.items():
    if rs in knee and rs in recon_knee:
        k=knee[rs]
        # already oriented to exposure allele via recon_knee
        dat_k.append(dict(beta_o=recon_knee[rs], se_o=k["se_ln"], beta_e=h["beta_e"], se_e=h["se_e"]))
bk, sek = ivw_fixed(dat_k)
bkeg, sekeg, ikeg, sikeg = mr_egger(dat_k)
print("\n=== Independent replication: BMD -> knee-OA (GCST90566800) ===")
print(f"n={len(dat_k)} IVW beta={bk:.5f} se={sek:.5f} OR/SD={math.exp(bk):.4f} [{math.exp(bk-1.96*sek):.4f},{math.exp(bk+1.96*sek):.4f}]")
print(f"MR-Egger slope={bkeg:.5f} se={sekeg:.5f} OR/SD={math.exp(bkeg):.4f} intercept={ikeg:.5f}")

# ---- save per-SNP table ----
per_snp=[]
for rs, h in inst.items():
    if rs in hip_oriented:
        d=hip_oriented[rs]
        per_snp.append(dict(rs=rs, beta_e=h["beta_e"], se_e=h["se_e"],
                            hip_lnOR=d["beta_o"], hip_se=d["se_o"], hip_OR=d["OR"],
                            hip_N=d["N"], knee_recon=recon_knee.get(rs), harm_beta_o=h["beta_o"]))
per_snp.sort(key=lambda x:x["rs"])

result = dict(
    validation=dict(n=len(validation), mean_diff=statistics.mean(diffs), max_abs_diff=max(abs(d) for d in diffs), pearson_r=r),
    F1=dict(n=len(dat), ivw_fixed=dict(beta=b_ivw, se=se_ivw, OR_sd=math.exp(b_ivw),
                ci_lo=math.exp(b_ivw-1.96*se_ivw), ci_hi=math.exp(b_ivw+1.96*se_ivw),
                p=2*(1-0.5*(1+math.erf(b_ivw/(se_ivw*math.sqrt(2)))))),
              ivw_random=dict(beta=b_ivw_r, se=se_ivw_r, OR_sd=math.exp(b_ivw_r), tau2=tau2),
              egger=dict(beta=b_eg, se=se_eg, OR_sd=math.exp(b_eg), intercept=int_eg, se_int=se_int),
              wmedian=dict(beta=b_wm, OR_sd=math.exp(b_wm))),
    knee_replication=dict(n=len(dat_k), ivw=dict(beta=bk, se=sek, OR_sd=math.exp(bk),
                ci_lo=math.exp(bk-1.96*sek), ci_hi=math.exp(bk+1.96*sek)),
                egger=dict(beta=bkeg, se=sekeg, OR_sd=math.exp(bkeg), intercept=ikeg, se_int=sikeg)),
    per_snp=per_snp,
)
json.dump(result, open(OUT,"w"), indent=2)
print("\nSaved", OUT)
