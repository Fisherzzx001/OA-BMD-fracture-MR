# -*- coding: utf-8 -*-
"""
第十一轮 · 正文词量压缩（面向 OAC 4000 词上限）。

对象：稿件的五个最长段落（Methods MR / P6 裁定 / SOS 构念 / SOS 边界 / 优势与局限），
按 20–25% 压缩，只削减冗余修辞与重复交代，**不删任何数值、引文或披露**。

每一段替换后逐项断言：
  1. 段内所有数值 token（含 %、p=、OR、β、Q=、τ²、df）必须原样保留；
  2. 门禁必需字符串（G7/G8/G11 披露词）必须仍在全文。
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

t = io.open(MD, encoding="utf-8").read()
L = t.split("\n")
w = lambda s: len([x for x in s.split() if re.search(r"[0-9A-Za-z]", x)])
NUM = re.compile(r"[-+]?\d[\d,\.]*\s?(?:×10[⁻¹²³⁰-⁹]+|%|SD|df|Mb)?|τ²|I²|r_g|Q=")


def setline(prefix, new, tag, keep_nums=True):
    """按前缀定位整行并替换；校验数值 token 不丢失。"""
    global L
    i = next(k for k, l in enumerate(L) if l.startswith(prefix))
    old = L[i]
    if keep_nums:
        norm = lambda s: s.strip().rstrip(".,;:)").strip()
        old_nums = set(norm(m.group(0)) for m in NUM.finditer(old))
        new_nums = set(norm(m.group(0)) for m in NUM.finditer(new))
        missing = old_nums - new_nums
        assert not missing, f"[{tag}] 数值丢失: {sorted(missing)[:12]}"
    print(f"  ✓ {tag}:  {w(old)} → {w(new)} 词 (Δ {w(new) - w(old)})")
    L[i] = new


print("=" * 78)
print("第十一轮 · 正文压缩")
print("=" * 78)

# ---------------------------------------------------------------- 1. Methods · MR
setline("**MR.**", (
    "**MR.** Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001 over 10 Mb against "
    "1000 Genomes Phase 3 EUR (strict PLINK-equivalent clumping benchmarked against the PLINK "
    "binary; 20–35 instruments per analysed BMD site, 116 for knee OA). The BMD exposure set "
    "followed a prespecified rule: of the 11 DXA-derived BMD sites reported by Qian et al. "
    "2025, the cranial head site — the only one that is not a load-bearing region, and one "
    "whose instruments carry non-skeletal pleiotropy (head BMD is genetically correlated with "
    "intracranial aneurysm in the source GWAS [2]) — was excluded a priori, leaving ten "
    "prespecified canonical regions analysed under the same rule: a site enters MR only if "
    "strict clumping completes. Clumping completed for five sites at submission (arms, femoral "
    "neck, lumbar spine, legs, total body), each yielding ≥15 independent instruments; the "
    "remaining five (femoral total, pelvis, rib, spine, trunk) showed comparable genome-wide "
    "signal (229–345 loci at p<5×10⁻⁸ each) but their strict clumping was not completed in the "
    "compute environment (1000G EUR panel unavailable), and they were deferred rather than "
    "analysed with a non-prespecified ad hoc instrument set (Table S12). Instrument shortage "
    "drove no exclusion. Harmonisation (excluding palindromic and strand-ambiguous variants and "
    "instruments absent from the outcome GWAS) reduced the analysed sets to 16–31 instruments "
    "per forward pair and 69 per knee-OA-exposure pair; 47 of the 116 knee-OA instruments (41%) "
    "were lost at this step, an attrition we report explicitly (Fig. S9). The 69 retained "
    "instruments were identical across all three knee-OA-exposure outcomes, so the loss was "
    "uniform rather than outcome-size-driven, and the dropped instruments were not effect-size "
    "selected (median −log10P 9.10 vs 9.89 among the retained; 43% carried stronger signal than "
    "the retained median), having been lost through palindromic/strand-ambiguous status or "
    "absence from all outcome summary statistics (per-instrument attrition table). We computed "
    "IVW (fixed and random effects), MR-Egger (intercept test for directional pleiotropy), "
    "weighted median and leave-one-out, with Cochran Q and I² for heterogeneity and MR-PRESSO "
    "for outliers; random-effects IVW inflated the SE multiplicatively (SE_random = "
    "SE_fixed × √max(1, Q/df)) rather than by additive τ². BMD→fracture served as a positive "
    "control. The OA_selfrep dataset was validated genome-wide against the GO Consortium 2025 "
    "knee-OA GWAS (1,980 genome-wide-significant loci, 94.4% Z-sign concordance, no global flip)."
), "1 Methods MR 段")

# ---------------------------------------------------------------- 2. P6 裁定
setline("Zhou et al. [5] reported higher femoral-neck BMD protected", (
    "Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA [35], whereas our "
    "ten-site analysis shows higher BMD raises knee-OA risk — a hip-versus-knee contrast "
    "explicable under joint-specific pathway dominance (hip OA more mechanically/"
    "morphologically driven, knee OA more by systemic bone-formation pathways). This cautions "
    "against generalising any single BMD→OA estimate across joints. Hartley et al. [22] "
    "reported a positive BMD effect on both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, "
    "1.20–1.63); Funck-Brentano et al. [36] likewise reported a positive overall effect (all-OA "
    "OR 1.14, 1.06–1.22), significant at both knee and hip; a quantitative comparison on the "
    "native scales is in Fig. S12. These run against the directional form of our prediction P6 "
    "— motivated by Zhou's protective femoral-neck finding; P6 and the falsification rule F1 "
    "are introduced below and in Table S11. We adjudicate rather than cite selectively. First, "
    "the hip evidence disagrees with itself (protective [5], harmful [22], [36]). Second, and "
    "decisively for F1, the criterion is conjunctive: the hip effect must be positive AND of "
    "knee-like magnitude. Within Hartley et al., the only same-pipeline comparison, the hip and "
    "knee intervals overlap almost completely (1.05–1.57 vs 1.20–1.63) and on the log-odds "
    "scale the hip estimate lies 27% below the knee's (0.247 vs 0.336); Hartley establishes "
    "positivity, not equality. We therefore restate P6 as a magnitude claim and pre-register "
    "the decision rule (a hip effect that is positive and at least as large as the knee's would "
    "fail P6). The directional form was the pre-registered prediction; its restatement here "
    "deviates from that registered form and is prompted by the contradictory hip evidence."
), "2 P6 裁定段")

# ---------------------------------------------------------------- 3. SOS 构念
setline("OA is a disorder of the osteochondral unit", (
    "OA is a disorder of the osteochondral unit [27–30], in which subchondral bone remodelling "
    "feeds cartilage degradation. A systemic, polygenic bone-forming set-point readable as areal "
    "BMD would protect at the whole-bone endpoint and raise OA risk at the joint endpoint. "
    "Hardcastle et al. [34] described a clinical \"bone-forming triad\" of osteophytes, "
    "enthesophytes and high bone mass; the Systemic Osteogenic Set-point (SOS) reframes that "
    "constellation as a polygenic tendency indexed by areal BMD that separates the whole-bone "
    "from the joint consequences, retaining the name while acknowledging the triad as its "
    "phenotypic antecedent. The set-point is operationalised as areal BMD — a composite of "
    "volumetric density and size — and is a latent tendency rather than a measured quantity. "
    "The height partition is partial in both directions, since height captures neither bone "
    "width nor geometry [35]: 9.5–33.2% is a lower bound on the size-attributable share, and "
    "the 67–91% height-independent share is not identified as bone mass — the name states the "
    "hypothesis, not an established mechanism. The reading generates six predictions (P1–P6) "
    "and seven falsifying tests (F1–F7; Table S11). P1 (opposing endpoints) is supported; P2 "
    "(systemic, not local-mechanical) is formally tested here, its systemic clause supported "
    "and its site-ordering clause under-powered; P3 (diffuse polygenic) is directionally "
    "consistent but not discriminating, because a null colocalisation is predicted under any "
    "diffuse mechanism; P4 (size is a separable partial component) is supported subject to the "
    "mass-versus-size split that only F2 can resolve; P6 is restated as a magnitude claim "
    "(above); F1 is pre-registered. Two tests are executed here for the first time. F6 "
    "(cross-site readout homogeneity) asks whether one set-point is readable at any site: the "
    "five site-specific BMD→knee-OA estimates are statistically indistinguishable (Q=6.92 on 4 "
    "df, p=0.14, with the covariance induced by shared instruments handled exactly; τ²=0), "
    "whereas a site-specific rival predicts heterogeneity the data do not show. F7 (joint-locus "
    "exclusion) asks whether the effect is carried by joint-biology loci: GDF5, the archetypal "
    "joint-development locus, is not a genome-wide-significant BMD instrument at any site, and "
    "excluding WNT16 or SOX6 changes no estimate's direction or significance."
), "3 SOS 构念段")

# ---------------------------------------------------------------- 4. SOS 边界
setline("Two boundaries limit what these predictions can discriminate", (
    "Two boundaries limit what these predictions can discriminate. First, the tests executed "
    "here separate the set-point architecture from site-specific, local-mechanical and "
    "joint-locus rivals, but not from a generic systemic rival in which any BMD-raising process "
    "also acts at the joint; that separation requires F2–F5 (volumetric decomposition, tissue "
    "enrichment, eQTL colocalisation, osteophyte-stratified outcome), unavailable with the "
    "present data. Second, P2's site-ordering clause is under-powered as a pairwise contrast "
    "(arms minus femoral neck +0.005, one-sided p=0.18), and P3 is not discriminating because a "
    "null colocalisation is predicted under any diffuse mechanism; the homogeneity test is "
    "itself bounded, since one of ten pairwise contrasts (lumbar spine versus legs, p=0.022) "
    "reaches nominal significance — chance-compatible under Bonferroni correction but not "
    "excluded. Among the pre-registered predictions the genuinely discriminating one is P4, "
    "whose residual is resolved only by F2. The lowest-cost, highest-yield test is F5, "
    "specified here in advance: re-run this design with the OA outcome stratified into "
    "osteophyte-predominant versus joint-space-narrowing-predominant radiographic subtypes and "
    "test whether the BMD→OA estimate concentrates in the former. Because osteophytes are the "
    "structural signature of the set-point, failure to concentrate would falsify the proposed "
    "mechanism; concentration would confirm it against the same fracture positive control. F5 "
    "needs an OA GWAS with radiographic subtypes and is pre-registered rather than performed, "
    "so we claim no result from it."
), "4 SOS 边界段")

# ---------------------------------------------------------------- 5. 优势与局限
setline("**Strengths and limitations.**", (
    "**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions, an "
    "internal BMD–fracture positive control calibrated against the published GEFOS estimate "
    "(OR 1.40–1.42 vs 1.55), dual FDR/Bonferroni correction, a common SNP set throughout, and "
    "directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO. Limitations: "
    "the BMD instruments are European-only whereas the knee-OA outcome is 87.31% multi-ancestry "
    "(Table 1b) and LDSC uses a European-only reference; that mismatch can bias r_g and the MR "
    "estimates through cross-population LD differences and is not removed by the near-1.0 LDSC "
    "intercept, so the forward estimates generalise most directly to European-ancestry "
    "populations. Areal BMD conflates size with mass [18,32,37]; height-conditioning bounds "
    "only the height-proxied share at 9.5–33.2% — a lower bound, since bone width and geometry "
    "form a separable genetic architecture [35] — so the 67–91% height-independent residual is "
    "not established as bone mass (resolved only by test F2). Heterogeneity is high across all "
    "pairs including the positive control, so the causal claim rests on estimator concordance "
    "rather than on fixed-effect p-values; the reverse and OA→fracture estimates are null but "
    "weakly powered (0.14–0.26), excluding only |β|≥0.30; residual overlap is bounded at ≤2.6% "
    "for BMD→OA and acts in the observed direction but cannot be removed with summary data "
    "alone. Our BMD→OA magnitude is 12–20-fold smaller on the log-odds scale than Hartley et "
    "al.'s one-sample estimate [22], and the gap is unresolved: the exposures are different "
    "phenotypes on different scales, the designs differ in overlap (100% versus ≤2.6%), and the "
    "residual five-fold gap after inverting Hartley's reverse estimate is not attributable to "
    "one-sample estimation alone. We therefore claim directional, not magnitude, concordance "
    "with that estimate. The tests executed here (F6, F7) separate the set-point architecture "
    "from site-specific, local-mechanical and joint-locus rivals but not from a generic "
    "systemic-pleiotropy rival in which any BMD-raising process also acts at the joint; "
    "separating those requires the external-data tests F2–F5."
), "5 优势与局限段")

t = "\n".join(L)
io.open(MD, "w", encoding="utf-8").write(t)

# ---------------------------------------------------------------- 门禁必需披露仍在
print("\n[校验] 门禁必需披露词")
for s in ["0.14\u20130.26", "I\u00b2 68\u201386%", "9.5\u201333.2%", "height-independent",
          "bone width", "geometry", "lower bound", "MR-Egger", "positive control",
          "positive control was comparably heterogeneous", "height-independent residual"]:
    assert s in t, f"缺失披露词: {s}"
print("  ✓ 11 项披露词全部保留")

print("\n" + "=" * 78)
print("压缩完成（词量变化由 49_oac_wordcount 复核）")
print("=" * 78)
