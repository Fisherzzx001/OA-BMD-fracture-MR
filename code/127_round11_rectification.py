# -*- coding: utf-8 -*-
"""
第十一轮整改：逐条落实第九轮预审报告 (docs/112) 的 P0/P1/P2 与报告内含的每一条批评。

本轮处理的条目（对照 112 报告）：
  Editor  #2  摘要以固定效应 p 领头，因果语气与正文 pleiotropy-robust 框架不一致 (L24)
  R3      #1  "first direct MR test of OA→any-fracture" 伪创新风险 (L92)
  R3      #2  [31] 误引（Tong 2025 踝足疾病→改为 Zhu 2025 关节置换位点特异性）(L32 + 参考文献)
  R3      #3  2023–2026 BMD→OA MR 文献补检（Zhu 2025 入库）
  R3      #4  Hartley 调和段压缩约 1/3 (L86)
  R2 建议 #4  5/10 位点选择逻辑：预设规则 + canonical DXA 位点 + "five analysed sites" 限定 (L48/L34/摘要/Highlights)
  R2 建议 #5  MR-Egger 功率限定 (L64)
  R2 建议 #6  femoral-neck→fracture Egger intercept p=0.018 正文点明 (L62)
  R4      #1  P6/F1 前向引用断裂，首次出现处加预告指针 (L84)
  R4      #3  方向一致 / 量级不一致的分层在摘要呼应
  格式门禁    Highlights 3–5 条 × ≤85 字符；摘要 ≤250 词
  陈旧值      Table 1 "N used" 列改标实际使用的 representative N=32,000 并补 ‡ 脚注

已在前序轮次完成、本轮仅复验不改：
  P0 MRlap 方向拆分 / P0 femoral-head(颅骨)排除 / P0 Highlights "explains" /
  P1 Effect-scale 断言式 / R1 SOS 鉴别性检验 F6+F7 / R1 forming-vs-size 严格术语
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

t = io.open(MD, encoding="utf-8").read()
orig = t
n_edits = 0


def rep(old, new, tag):
    """精确替换一次，带断言。"""
    global t, n_edits
    n = t.count(old)
    assert n == 1, f"[{tag}] 期望 1 处，实际 {n} 处 :: {old[:90]}"
    t = t.replace(old, new)
    n_edits += 1
    print(f"  ✓ {tag}")


def repline(idx, new, tag, must_startswith=None, must_endswith=None):
    """整行替换（0-based 行号），带首尾断言。"""
    global t, n_edits
    L = t.split("\n")
    old = L[idx]
    if must_startswith:
        assert old.startswith(must_startswith), f"[{tag}] 行首不符: {old[:80]}"
    if must_endswith:
        assert old.endswith(must_endswith), f"[{tag}] 行尾不符: ...{old[-60:]}"
    L[idx] = new
    t = "\n".join(L)
    n_edits += 1
    print(f"  ✓ {tag}")


print("=" * 78)
print("第十一轮整改 · 稿件补丁")
print("=" * 78)

# ---------------------------------------------------------------- 1. 摘要
print("\n[1] 摘要：因果语气对齐 + 方向/量级分层呼应 + 词量压缩")

rep("using linkage-disequilibrium score regression across 78 phenotype pairs and "
    "two-sample Mendelian randomisation (IVW, MR-Egger, weighted median, MR-PRESSO, "
    "leave-one-out), with BMD's protection against fracture as a positive control.",
    "using linkage-disequilibrium score regression across 78 phenotype pairs and "
    "two-sample Mendelian randomisation, with BMD→fracture as a positive control.",
    "1a 摘要方法句去掉与结果重复的估计量清单")

rep("both OA definitions correlated positively with every BMD site (r_g 0.12–0.25) "
    "and with each other (r_g 0.84).",
    "both OA definitions correlated positively with every BMD site (r_g 0.12–0.25).",
    "1b 摘要去掉次要结果 r_g 0.84")

rep("and showed genetically higher BMD raising knee-OA risk (β +0.017 to +0.028 per SD; "
    "p_IVW ≤6.4×10⁻⁴). Heterogeneity was high (I² 68–86%) but matched the positive "
    "control (78–83%), and all five estimates kept a positive direction across IVW, "
    "MR-Egger (intercepts p>0.24), weighted median and MR-PRESSO; 2 of 5 remained "
    "significant under random effects (3 of 5 under DerSimonian–Laird).",
    "and showed genetically higher BMD raising knee-OA risk (β +0.017 to +0.028 per SD; "
    "fixed-effect p ≤6.4×10⁻⁴), a direction preserved across MR-Egger (intercepts "
    "p>0.24), weighted median and MR-PRESSO. Heterogeneity was high (I² 68–86%) but "
    "matched the positive control (78–83%); 2 of 5 estimates remained significant under "
    "random effects (3 of 5 under DerSimonian–Laird), so the causal reading rests on "
    "estimator concordance rather than on fixed-effect p-values.",
    "1c 摘要因果语气：以估计量一致领头，固定效应 p 降为附证")

rep("Conditioning on height removed 9.5–33.2% of the effect, leaving 67–91% "
    "height-independent, a residual not resolved into bone mass versus geometry.",
    "Conditioning on height removed 9.5–33.2%; the 67–91% height-independent residual "
    "is not resolved into bone mass versus geometry.",
    "1d 摘要 size 句紧缩")

rep("Reverse MR was null and underpowered (power 0.14–0.26).",
    "Reverse estimates were non-significant and underpowered (power 0.14–0.26).",
    "1e 摘要反向句紧缩")

rep("One BMD-raising profile thus protects against fracture yet raises knee-OA risk; "
    "the asymmetry rests on the forward estimates, whose overlap-induced bias is bounded "
    "at ≤2.6% for BMD→OA and 12% for BMD→fracture and acts in the direction of r_g;",
    "One BMD-raising profile thus protects against fracture yet raises knee-OA risk — "
    "directionally concordant with, but 12–20-fold smaller on the log-odds scale than, "
    "the largest published one-sample estimate — and the asymmetry rests on forward "
    "estimates whose overlap-induced bias is bounded at ≤2.6% (BMD→OA) and 12% "
    "(BMD→fracture), acting in the direction of r_g;",
    "1f 摘要呼应方向一致/量级不一致分层")

# ---------------------------------------------------------------- 2. Highlights 格式
print("\n[2] Highlights：3–5 条 × ≤85 字符（OAC 硬约束）")

_NEW_HL = (
    "- Higher genetically predicted BMD raises knee OA risk yet protects against fracture.\n"
    "- BMD→OA was positive at all five analysed sites; 2 of 5 held under random effects.\n"
    "- 9.5–33.2% of the BMD→OA effect is attributable to the height-proxied size pathway.\n"
    "- Reverse estimates were non-significant and underpowered (power 0.14–0.26).\n"
    "- A systemic bone-forming set-point remains a hypothesis consistent with these data."
)
L = t.split("\n")
_h0 = next(i for i, l in enumerate(L) if l.strip() == "## Highlights")
_h1 = next(i for i, l in enumerate(L) if l.startswith("**Keywords"))
_hl_idx = [i for i in range(_h0, _h1) if L[i].startswith("- ")]
assert len(_hl_idx) == 6, f"Highlights 应为 6 条旧条目，实际 {len(_hl_idx)}"
assert _hl_idx == list(range(_hl_idx[0], _hl_idx[0] + 6)), "Highlights 应连续"
L[_hl_idx[0]:_hl_idx[-1] + 1] = _NEW_HL.split("\n")
t = "\n".join(L)
n_edits += 1
print("  ✓ 2 Highlights 重写为 5 条")
for _h in _NEW_HL.split("\n"):
    _c = len(_h) - 2
    print(f"      {_c:>3} chars {'OK' if _c <= 85 else 'OVER'}  {_h[2:]}")

# ---------------------------------------------------------------- 3. R3#2 误引
print("\n[3] R3#2：[31] 误引替换为位点特异性证据 (Zhu 2025)")

rep("other joint-specific MR analyses of OA report a similar cautionary pattern [31].",
    "and a 2025 two-sample MR of BMD against OA progression reports the same site-by-joint "
    "dependence, with femoral-neck BMD raising hip-replacement risk (OR 1.27, 95% CI "
    "1.13–1.43) but not knee-replacement risk [31].",
    "3a L32 引用替换为关节特异性实证")

_NEW31 = ("31. Zhu R, Xing X, Bian J, Zhang X, Ge L, Cai G. Causal association between bone "
          "mineral density and the risk of joint replacement in patients with osteoarthritis: "
          "a Mendelian randomization study. *Clin Rheumatol*. 2025;44(2):823–830. "
          "doi:10.1007/s10067-024-07289-5. PMID 39738847.")
t2, cnt = re.subn(r"^31\. .*$", lambda m: _NEW31, t, count=1, flags=re.M)
assert cnt == 1, f"[3b] 参考文献 31 替换失败 ({cnt})"
t = t2
n_edits += 1
print("  ✓ 3b 参考文献 [31] 换为 Zhu et al. 2025 (PMID 39738847)")
assert "Tong X" not in t, "Tong 2025 应已移除"

# ---------------------------------------------------------------- 4. R2#4 位点选择
print("\n[4] R2#4：5/10 位点选择的预设性说明 + analysed-sites 限定")

rep("Of these, strict PLINK-equivalent clumping was completed for five sites at submission "
    "(arms, femoral neck, lumbar spine, legs, total body), each yielding ≥15 independent "
    "instruments and entering MR;",
    "All ten are prespecified canonical DXA regions analysed under one prespecified rule: a "
    "site enters MR only if strict clumping completes at p<5×10⁻⁸, r²<0.001 over 10 Mb. "
    "Strict clumping was completed for five sites at submission (arms, femoral neck, lumbar "
    "spine, legs, total body), each yielding ≥15 independent instruments and entering MR;",
    "4a L48 预设规则 + canonical DXA 位点")

rep("(i) ten regional BMD sites rather than a single composite BMD, so site specificity is "
    "tested directly;",
    "(i) ten regional BMD sites rather than a single composite BMD, so site specificity is "
    "tested directly at the sites that entered MR (five of ten; §2.3);",
    "4b L34 位点特异性限定到进入 MR 的五个位点")

# ---------------------------------------------------------------- 5. R2#6 阳性对照
print("\n[5] R2#6：femoral-neck→fracture 的 Egger intercept p=0.018 正文点明")

rep("The femoral-neck set (25 SNPs) gave a concordant result (β=−0.351; OR 1.42, 1.34–1.51); "
    "the conservative 1-Mb set gave β=−0.366 (OR 1.44).",
    "The femoral-neck set (25 SNPs) gave a concordant result (β=−0.351; OR 1.42, 1.34–1.51), "
    "but its MR-Egger intercept was not null (p=0.018), a directional-pleiotropy signal we "
    "report rather than absorb: its pleiotropy-robust slope (−0.343) and weighted median "
    "(−0.365) nevertheless agree with IVW, so the control calibrates the pipeline while its "
    "own diagnostic is not perfectly clean. The conservative 1-Mb set gave β=−0.366 "
    "(OR 1.44).",
    "5 L62 阳性对照的异质性信号")

# ---------------------------------------------------------------- 6. R2#5 Egger 功率
print("\n[6] R2#5：MR-Egger intercept 的功率限定")

rep("MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign.",
    "MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign; with "
    "16–31 instruments the intercept test has limited power, so a non-significant intercept "
    "excludes unbalanced pleiotropy of detectable size rather than pleiotropy itself.",
    "6 L64 Egger 功率限定")

# ---------------------------------------------------------------- 7. R4#1 前向引用
print("\n[7] R4#1：P6/F1 前向引用断裂")

rep("These run against the directional form of our prediction P6 (motivated by Zhou's "
    "protective femoral-neck finding).",
    "These run against the directional form of our prediction P6 — motivated by Zhou's "
    "protective femoral-neck finding; P6 and the falsification rule F1 are introduced below "
    "and in Table S11.",
    "7 L84 P6/F1 首次出现处加预告指针")

# ---------------------------------------------------------------- 8. R3#1 首例声明
print("\n[8] R3#1：'first' 创新声明软化")

rep("This is the first direct MR test of OA→any-fracture.",
    "To our knowledge this is among the first direct MR tests of OA→any-fracture.",
    "8 L92 首例声明软化")

# ---------------------------------------------------------------- 9. R3#4 Hartley 压缩
print("\n[9] R3#4：Hartley 调和段压缩约 1/3")

_NEW86 = (
    "The directional agreement with Hartley et al. [22] is matched by a magnitude "
    "disagreement we state rather than leave to be discovered: their knee estimate is OR 1.40 "
    "(1.20–1.63) per SD higher eBMD — 0.336 on the log-odds scale — against 0.017–0.028 "
    "here, a 12–20-fold difference. Three features make this a difference of scale and "
    "design rather than of measured biology. First, the exposures differ: Hartley "
    "instrumented heel-ultrasound eBMD, a raw attenuation measure that scales with bone size, "
    "whereas we instrumented DXA areal BMD, mineral content divided by projected area and "
    "therefore area-normalised, so the two per-SD effects are not denominated in a common "
    "unit and their ratio is not a disagreement about the causal effect per unit of bone. "
    "Second, the designs differ: Hartley's headline estimate is one-sample two-stage least "
    "squares on unweighted allele scores measured in the same individuals as the outcome, so "
    "overlap is 100% by construction, whereas ours is bounded at ≤2.6% and an MRlap "
    "correction for overlap, weak-instrument and winner's-curse bias moved every estimate a "
    "median 20.5% further from the null. Third, our pipeline is calibrated rather than "
    "asserted: the BMD→fracture positive control returns OR 1.40–1.42 against the published "
    "GEFOS 1.55, an 8.4–9.6% shortfall on the ratio scale, so a systematic deflation would "
    "have to be specific to the OA outcome. Inverting Hartley's reverse estimate onto the "
    "forward scale (β 4.16 SD per doubling of knee-OA risk, implying OR ≈1.18 per SD higher "
    "eBMD) leaves a residual gap of roughly five-fold, so the discrepancy is not an artefact "
    "of one-sample estimation alone. We do not claim that our estimate supersedes theirs: the "
    "two studies instrument different phenotypes under different designs, and no analysis "
    "available to us places them on a common scale. We read the directions as concordant and "
    "the magnitudes as unresolved, and record the disagreement as a limitation."
)
_L = t.split("\n")
_ix86 = next(i for i, l in enumerate(_L)
             if l.startswith("The directional agreement with Hartley"))
repline(_ix86, _NEW86, "9 L86 Hartley 段压缩",
        must_startswith="The directional agreement with Hartley",
        must_endswith="record the disagreement as a limitation.")
_old_words = len([w for w in orig.split("\n")[_ix86 + 1].split()
                  if re.search(r"[0-9A-Za-z]", w)])
_new_words = len([w for w in _NEW86.split() if re.search(r"[0-9A-Za-z]", w)])
print(f"      L86 词数 {_old_words} → {_new_words} (Δ {_new_words - _old_words})")

# ---------------------------------------------------------------- 10. 陈旧值：Table 1 N
print("\n[10] 陈旧值清理：Table 1 'N used' 改标实际使用的 representative N=32,000")

L = t.split("\n")
_nfix = 0
for i, l in enumerate(L):
    if l.startswith("| BMD") and "31,873‡" in l:
        L[i] = l.replace("31,873‡", "32,000‡")
        _nfix += 1
assert _nfix == 10, f"Table 1 应有 10 行 BMD，实际替换 {_nfix} 行"
t = "\n".join(L)
n_edits += 1
print(f"  ✓ 10a Table 1 BMD 行 N used 31,873‡ → 32,000‡（{_nfix} 行）")

rep("§ OA_selfrep N confirmed at 462,933 (38,472 cases / 424,461 controls; IEU OpenGWAS "
    "`ukb-b-14486`).",
    "§ OA_selfrep N confirmed at 462,933 (38,472 cases / 424,461 controls; IEU OpenGWAS "
    "`ukb-b-14486`). ‡ BMD h² was computed at a rounded representative N of 32,000 (the "
    "per-site GWAS N is 31,873); because h² ∝ 1/N the rounding makes the reported h² 0.40% "
    "lower than at the exact N, an effect far smaller than the reported SEs.",
    "10b 补 ‡ 脚注说明 representative N")

# ---------------------------------------------------------------- 写回 + 自检
assert t != orig, "文件未发生变化"
io.open(MD, "w", encoding="utf-8").write(t)

print("\n" + "=" * 78)
print(f"稿件补丁完成：{n_edits} 处修改")

# 自检：不得残留的表述
print("\n[自检] 整改后残留检查")
_bad = []
for _w in ["This is the first direct MR test", "Tong X", "p_IVW ≤6.4×10⁻⁴",
           "other joint-specific MR analyses"]:
    if _w in t:
        _bad.append(_w)
assert not _bad, f"仍残留旧表述: {_bad}"
print("  ✓ 旧表述零残留")

for _w in ["To our knowledge this is among the first", "p=0.018", "limited power",
           "prespecified canonical DXA regions", "introduced below", "32,000",
           "12–20-fold smaller", "fixed-effect p ≤6.4×10⁻⁴"]:
    assert _w in t, f"新表述缺失: {_w}"
print("  ✓ 新表述全部落盘")

_hl = t.split("## Highlights")[1].split("**Keywords")[0]
_hls = [l[2:].strip() for l in _hl.split("\n") if l.startswith("- ")]
assert 3 <= len(_hls) <= 5, f"Highlights 条数 {len(_hls)} 不在 3–5"
for h in _hls:
    assert len(h) <= 85, f"Highlights 超 85 字符 ({len(h)}): {h}"
print(f"  ✓ Highlights {len(_hls)} 条，最长 {max(len(h) for h in _hls)} 字符")

_abs = t.split("## Abstract")[1].split("##")[0]
_aw = len([w for w in _abs.split() if re.search(r"[0-9A-Za-z]", w)])
print(f"  ✓ 摘要词数 {_aw}（OAC 上限 250）{'OK' if _aw <= 250 else 'OVER'}")

_body = 0
for _s in re.split(r"\n#{1,6}\s", t):
    pass
print("=" * 78)
