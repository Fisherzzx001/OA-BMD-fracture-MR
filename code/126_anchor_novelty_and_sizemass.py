# -*- coding: utf-8 -*-
"""
126 · 第十轮深挖落地：(1) opposing-endpoints 新意锚定 (2) bone-forming 与 size/geometry 未分离

对应第九轮预审报告 R1 致命漏洞 2/3：
  漏洞3 "opposing endpoints 新意被高估" -> 贡献锚定到"单一遗传工具集 + 单一校准管线的联合证明"，
       而非对立概念本身（概念自 Foss&Byers 1972 / Dequeker 2003 即已存在）。
  漏洞2 "bone-forming 与 bone-size/shape set-point 未分离" -> 67-91% 严格表述为 height-independent，
       明确不等于 mass；height 不含骨宽/几何（Faber [35] 可分离遗传架构），9.5-33.2% 是 size 归因下界；
       残差的 mass-vs-geometry 鉴别只能由 F2 完成（本地无 OA/fracture 全基因组 sumstats，无法补做）。

所有替换带断言；不改任何已报告数值。
"""
import io, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

text = io.open(MD, encoding="utf-8").read()
orig_words = len([t for t in text.split() if re.search(r"[0-9A-Za-z]", t)])

EDITS = [
    # ---- Highlights H4：size 归因严格化（去掉 "only ... explained by skeletal size" 的 mass 暗示）----
    ("- BMD\u2192OA effect is only 9.5\u201333.2% explained by skeletal size (height proxy).",
     "- 9.5\u201333.2% of the BMD\u2192OA effect is attributable to the height-proxied size pathway."),

    # ---- Highlights H6：explains -> hypothesis（第九轮 P0-3 遗留项，同句落地）----
    ("- A systemic bone-forming set-point explains the opposing endpoint consequences.",
     "- A systemic bone-forming set-point remains a hypothesis consistent with these data."),

    # ---- 摘要：size 句严格化；MRlap 句压缩抵偿词量 ----
    ("Conditioning on height removed 9.5\u201333.2% of the effect, leaving 67\u201391% independent of the height\u2013size pathway.",
     "Conditioning on height removed 9.5\u201333.2% of the effect, leaving 67\u201391% independent of height, a residual not resolved into bone mass versus skeletal geometry."),
    ("an MRlap correction addressing overlap together with weak-instrument and winner's-curse bias moved all ten estimates further from the null (median 20.5%) without reversing any direction or conclusion.",
     "an MRlap correction for overlap, weak-instrument and winner's-curse bias moved all ten estimates further from the null (median 20.5%) without reversing any conclusion."),

    # ---- Intro gap 3：承认对立概念古老，锚定联合估计为新意 ----
    ("Third, the OA\u2013BMD\u2013fracture triangle has not been analysed jointly.",
     "Third, although each association is individually familiar [11,12,33], the OA\u2013BMD\u2013fracture triangle has not been analysed jointly."),

    # ---- Results L76：67-91% -> height-independent；残差不等于 mass ----
    ("leaving 67\u201391% independent of the height\u2013size pathway. The positive control confirmed a true mass effect survives conditioning: BMD\u2192fracture moved only \u22121.7% to \u22123.0%.",
     "leaving 67\u201391% independent of height. Because height indexes stature but not bone width or geometry [35], this residual is not identified as bone mass; the mass-versus-geometry split is resolved only by test F2. The positive control confirmed that the established BMD\u2192fracture effect survives the same conditioning: BMD\u2192fracture moved only \u22121.7% to \u22123.0%."),

    # ---- Discussion L82：全段重写（新意锚定 + forming/size 未分离披露 + Faber 定性支撑）----
    ("The paradox resolves once the two endpoints are separated. High BMD protects against fracture, confirmed here as a positive control, yet the same variants that raise BMD also raise OA risk. Subchondral bone sclerosis and osteophyte formation are hallmarks of OA and may share regulatory programs with systemic bone accrual [27\u201330]. That BMD protects against fracture yet raises OA risk implies a single BMD-increasing profile has opposing consequences at the two endpoints. One non-exclusive explanation is that the BMD\u2192OA effect is partly driven by skeletal size rather than bone mass, and that contribution is now quantified: areal BMD (DXA) is bone mineral content divided by projected area and is confounded by skeletal size [18,32]. MVMR conditioning on height attenuated the effect by 9.5\u201333.2%, while the fracture positive control moved only \u22121.7% to \u22123.0%, and height itself carried a direct OA-raising effect consistent with established epidemiology [24,25,26]. The height-independent share (67\u201391%) is therefore an upper bound on the mass-attributable component, pending a direct volumetric decomposition (test F2).",
     "The opposing-endpoint pattern is an old clinical observation [11,12]; what is new is its demonstration for a single genetic instrument set within one calibrated pipeline. One non-exclusive explanation of the OA signal is partial mediation by skeletal size rather than bone mass [18,32]. Conditioning on height attenuated the effect by 9.5\u201333.2% (the fracture control moved only \u22121.7% to \u22123.0%), and height carried a direct OA-raising effect consistent with established epidemiology [24,25,26]. Because height indexes stature but not bone width, cortical geometry or shape [35], the 9.5\u201333.2% is a lower bound on the size-attributable share, and the height-independent residual (67\u201391%) is not identified as bone mass; its mass-versus-geometry split is resolved only by test F2. Consistent with a mass component without establishing it, the fracture control is height-transparent, and hip shape shows no causal effect on hip OA in the largest skeletal-geometry GWAS [35]."),

    # ---- Discussion L88：Hardcastle/SOS 句压缩；operationalisation 句补 forming/size 边界 ----
    ("Hardcastle et al. [34] described a clinical \"bone-forming triad\" of osteophytes, enthesophytes and high bone mass in OA; the Systemic Osteogenic Set-point (SOS) reframes that local radiographic constellation as a polygenic, systemic tendency indexed by areal BMD that separates the whole-bone (fracture-protective) from the joint (OA-risk) consequences, extending the triad to a joint-versus-whole-bone split. We retain the SOS nomenclature for the construct while acknowledging the triad as its phenotypic antecedent.",
     "Hardcastle et al. [34] described a clinical \"bone-forming triad\" of osteophytes, enthesophytes and high bone mass in OA; the Systemic Osteogenic Set-point (SOS) reframes that constellation as a polygenic, systemic tendency indexed by areal BMD that separates the whole-bone (fracture-protective) from the joint (OA-risk) consequences, and we retain the name while acknowledging the triad as its phenotypic antecedent."),
    ("The set-point is operationalised as areal BMD, a composite of volumetric density and size that we partition only partially here (67\u201391% height-independent). It is a latent tendency indexed by areal BMD, not a directly measured quantity.",
     "The set-point is operationalised as areal BMD, a composite of volumetric density and size, and is a latent tendency rather than a measured quantity. The height partition is partial in both directions, since height captures neither bone width nor geometry [35]: 9.5\u201333.2% is a lower bound on the size-attributable share, and the 67\u201391% height-independent share is not identified as bone mass \u2014 the name states the hypothesis, not an established mechanism."),

    # ---- Limitations L94：size 下界 + 残差非 mass + 几何分解不可得 ----
    ("areal BMD conflates size with mass [18,32,37], and although height-conditioning bounds the size-attributable share at 9.5\u201333.2%, a direct volumetric decomposition (test F2) is not yet available.",
     "areal BMD conflates size with mass [18,32,37]; height-conditioning bounds only the height-proxied share of size at 9.5\u201333.2% \u2014 a lower bound, since bone width and geometry form a separable genetic architecture [35] \u2014 so the 67\u201391% height-independent residual is not established as bone mass (volumetric or geometric decomposition, test F2, unavailable)."),

    # ---- Conclusion L96：independent of skeletal size -> independent of height ----
    ("The effect pertains to the clinically measured areal-BMD phenotype, 67\u201391% of it independent of skeletal size; reverse causation",
     "The effect pertains to the clinically measured areal-BMD phenotype, 67\u201391% of it independent of height (not resolved into bone mass versus skeletal geometry); reverse causation"),

    # ---- S10 图题（figure legends，不计入 OAC 词数）----
    ("**S10.** Multivariable MR decomposition of the BMD\u2192OA effect conditional on height: (a) univariable vs MVMR estimates; (b) attenuation percentages with the BMD\u2192fracture positive control.",
     "**S10.** Multivariable MR decomposition of the BMD\u2192OA effect conditional on height: (a) univariable vs MVMR estimates; (b) attenuation percentages with the BMD\u2192fracture positive control; the height-independent residual is not resolved into bone mass versus skeletal geometry."),

    # ---- Table S9 脚注（表格，不计入）----
    ("Size-attributable fraction 9.5\u201333.2% across the two tabulated baselines (13.5\u201333.2% vs Table 3 and 9.5\u201321.2% vs the MVMR-set baseline); 67\u201391% independent of the height\u2013size pathway. Positive control moved \u22121.7% to \u22123.0% (mass effect intact).",
     "Height-attributable fraction 9.5\u201333.2% across the two tabulated baselines (13.5\u201333.2% vs Table 3 and 9.5\u201321.2% vs the MVMR-set baseline), a lower bound for total size; the 67\u201391% height-independent residual is not resolved into mass versus geometry. Positive control moved \u22121.7% to \u22123.0% (control intact)."),
]

fails = []
for i, (old, new) in enumerate(EDITS, 1):
    n = text.count(old)
    if n != 1:
        fails.append(f"R{i}: found {n} occurrences (expected 1)\n  OLD: {old[:110]}")
        continue
    text = text.replace(old, new)
    print(f"  R{i:2d} OK  ({len(old.split())}w -> {len(new.split())}w)")

if fails:
    print("\nFAILED:")
    for f in fails:
        print(" ", f)
    sys.exit(1)

io.open(MD, "w", encoding="utf-8").write(text)
new_words = len([t for t in text.split() if re.search(r"[0-9A-Za-z]", t)])
print(f"\nTotal document words: {orig_words} -> {new_words} ({new_words - orig_words:+d})")

# ---- 过度断言残留自检 ----
print("\n[over-claim residue scan]")
for bad in ["mass-attributable", "explains the opposing", "independent of skeletal size",
            "explained by skeletal size", "height\u2013size pathway", "mass effect intact",
            "true mass effect survives"]:
    hits = [m.start() for m in re.finditer(re.escape(bad), text)]
    status = "CLEAN" if not hits else f"RESIDUE x{len(hits)}"
    print(f"  {status:12s} {bad!r}")
    if hits:
        sys.exit(2)

# ---- 必备披露自检 ----
print("\n[required disclosures]")
for must in ["height-independent", "bone width", "geometry", "lower bound",
             "not identified as bone mass", "test F2", "[11,12]"]:
    ok = must in text
    print(f"  {'OK ' if ok else 'MISSING':8s} {must!r}")
    if not ok:
        sys.exit(3)
print("\nPATCH 126 COMPLETE")
