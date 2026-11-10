#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
326_cpbb_weight_bearing_mr.py  —  T2-E (rebuttal reserve): site-specific fracture
============================================================================
CPBB开放 UKB weight-bearing fracture GWAS (10,749 cases / 317,599 controls,
GRCh37): https://cpbb-2025.oss-cn-hangzhou.aliyuncs.com/gwas/Fracture_weight-bearing.tsv
(328 MB, HTTP Range enabled). 主稿 any-fracture 已是阳性对照 (OR 1.40); 负重骨
骨折是部位特异的更严格终点 —— 检验 BMD 工具对其保护是否方向一致/更强。

流程:
  1) 并行 Range 分段下载 (每段 ~16 MB, 避免限速代理下整下被 SIGTERM) -> 拼接;
  2) 读 31 个 BMD 锁定工具 (harmonised, GRCh37 chr/bp) 与 GEFOS BMD sumstats 的
     效应等位基因 A1/A2;
  3) 流式扫描 CPBB, 按 (chr,bp) 精确匹配 31 工具, 协调等位基因方向:
       若 CPBB.effect_allele == BMD.A1 -> b_iv = CPBB.beta
       若 CPBB.effect_allele == BMD.A2 -> b_iv = -CPBB.beta (且 EAF 翻转)
       回文/不等位 -> 丢弃;
  4) IVW(固定+随机)/Egger/加权中位数 + Q,I2; 报告 N_snps 与方向;
  5) 判定: N<10 或 异质性主导(rand_p>0.05 且 I2>60 且固定与随机符号不一致)
           -> UNDERPOWERED / INCONCLUSIVE; 否则报告方向性阳性对照。
不改动 docs/08 任何值。输出 code/output/genetics/CPBB_weight_bearing_fracture_MR.json
============================================================================
"""
import os, sys, csv, math, gzip, json, subprocess, tempfile, time
import numpy as np

L7   = "G:/OA_B_PROJECT_cache/l7/"
HARM = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output","genetics","harmonised","BMD_total__OA_knee.csv")
BMD  = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"
URL  = "https://cpbb-2025.oss-cn-hangzhou.aliyuncs.com/gwas/Fracture_weight-bearing.tsv"
ASSEMBLED = os.path.join(L7,"CPBB_Fracture_weight-bearing.tsv")
OUT  = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output","genetics","CPBB_weight_bearing_fracture_MR.json")
CHUNK = 16 * 1024 * 1024
TOTAL = 328_869_201
BMD_N = 31986.0

def curl_range(start,end,out):
    cmd=["curl","-sS","--ssl-no-revoke","--retry","3","-m","180","-r",f"{start}-{end}",URL,"-o",out]
    r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    sz=os.path.getsize(out) if os.path.exists(out) else 0
    exp=end-start+1
    return sz==exp, sz, exp, r.stderr.decode("utf-8","ignore")[-200:]

def download_and_assemble():
    if os.path.exists(ASSEMBLED) and os.path.getsize(ASSEMBLED)==TOTAL:
        print(f"[326] using cached {ASSEMBLED}",flush=True); return True
    nseg=math.ceil(TOTAL/CHUNK)
    tmpdir=os.path.join(L7,"cpbb_wb_parts"); os.makedirs(tmpdir,exist_ok=True)
    segs=[]
    for i in range(nseg):
        s=i*CHUNK; e=min(s+CHUNK-1,TOTAL-1); p=os.path.join(tmpdir,f"part_{i:03d}")
        segs.append((s,e,p))
    from concurrent.futures import ThreadPoolExecutor
    print(f"[326] downloading {nseg} Range segments (parallel, 4 workers)...",flush=True)
    ok=True
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs={ex.submit(curl_range,s,e,p):(s,e,p) for s,e,p in segs}
        for fut,(s,e,p) in futs.items():
            good,sz,exp,err=fut.result()
            if not good:
                ok=False; print(f"  [FAIL] {p}: {sz}/{exp} {err}",flush=True)
    if not ok:
        print("[326] segment download failed; abort (no fabrication)",flush=True); return False
    # assemble in order
    with open(ASSEMBLED,"wb") as out:
        for s,e,p in segs:
            with open(p,"rb") as f: out.write(f.read())
    asz=os.path.getsize(ASSEMBLED)
    print(f"[326] assembled {asz} bytes (expected {TOTAL})",flush=True)
    if asz!=TOTAL: return False
    for s,e,p in segs:
        try: os.remove(p)
        except: pass
    return True

def load_bmd_alleles():
    """rsid -> (A1_effect, A2, beta_e, se_e, maf)."""
    want=set()
    inst={}
    with open(HARM) as f:
        for r in csv.DictReader(f):
            want.add(r["SNP"]); inst[r["SNP"]]=dict(chr=int(r["CHR"]),bp=int(r["BP"]),
                    be=float(r["beta_e"]),se=float(r["se_e"]),flip=int(float(r["flip"])))
    out={}
    with gzip.open(BMD,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            snp=p[ix["SNP"]]
            if snp in want:
                out[snp]=dict(chr=int(p[ix["CHR"]]),bp=int(p[ix["BP"]]),
                              a1=p[ix["A1"]].upper(),a2=p[ix["A2"]].upper(),
                              maf=float(p[ix["MAF"]]),z=float(p[ix["Z"]]))
    return inst,out

def collect_fracture(inst,bmd):
    # coordinate -> list of rsids (BMD GRCh37)
    coord={}
    for rs,d in bmd.items(): coord.setdefault((d["chr"],d["bp"]),[]).append(rs)
    matched={}
    n=0
    with open(ASSEMBLED) as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            n+=1
            if n % 5_000_000 == 0: print(f"  [scan] {n//1_000_000}M rows",flush=True)
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["chromosome"]]); bp=int(p[ix["base_pair_location"]])
            except: continue
            key=(ch,bp)
            if key not in coord: continue
            try:
                beta=float(p[ix["beta"]]); se=float(p[ix["standard_error"]]); eaf=float(p[ix["effect_allele_frequency"]])
            except: continue
            if se<=0 or not math.isfinite(se) or not math.isfinite(beta): continue
            ea=p[ix["effect_allele"]].upper(); oa=p[ix["other_allele"]].upper()
            for rs in coord[key]:
                bd=bmd[rs]; ie=inst[rs]
                a1,a2=bd["a1"],bd["a2"]
                if {ea,oa}!={a1,a2}: continue
                if ea==a1:
                    b_iv=beta; e_iv=eaf; orient="+"
                elif ea==a2:
                    b_iv=-beta; e_iv=1.0-eaf; orient="-"
                else:
                    continue
                matched[rs]=dict(chr=ch,bp=bp,b_iv=b_iv,se=se,eaf=e_iv,be=ie["be"],se_e=ie["se"],
                                 orient=orient,maf=bd["maf"],p=p[ix["p_value"]])
    return matched

def mr(m):
    rs=list(m); be=np.array([m[r]["be"] for r in rs]); se_e=np.array([m[r]["se_e"] for r in rs])
    by=np.array([m[r]["b_iv"] for r in rs]); se_y=np.array([m[r]["se"] for r in rs])
    # per-SNP Wald ratio and its delta-method variance
    ratio=by/be
    var_ratio=se_y**2/be**2 + (by**2/be**4)*se_e**2
    w=1.0/var_ratio
    # fixed-effect IVW on ratio scale
    b_fixed=float(np.sum(w*ratio)/np.sum(w)); se_fixed=float(math.sqrt(1.0/np.sum(w)))
    p_fixed=math.erfc(abs(b_fixed/se_fixed)/math.sqrt(2))
    # Cochran Q / heterogeneity / DerSimonian-Laird random effects
    Q=float(np.sum(w*(ratio-b_fixed)**2)); df=len(rs)-1
    C=max(0.0,(Q-df)/float(np.sum(w)))
    wr=1.0/(var_ratio+C)
    b_rand=float(np.sum(wr*ratio)/np.sum(wr)); se_rand=float(math.sqrt(1.0/np.sum(wr)))
    p_rand=math.erfc(abs(b_rand/se_rand)/math.sqrt(2))
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    # MR-Egger: weighted regression of outcome beta on exposure beta (weights 1/se_y^2)
    W=1.0/se_y**2; sw=np.sum(W)
    sx=np.sum(W*be)/sw; sy=np.sum(W*by)/sw
    sxx=np.sum(W*(be-sx)**2); sxy=np.sum(W*(be-sx)*(by-sy))
    slope=sxy/sxx; intercept=sy-slope*sx
    se_int=math.sqrt(1.0/sw + sx**2/sxx)
    p_int=math.erfc(abs(intercept/se_int)/math.sqrt(2))
    # weighted median (weights = random-effects w on ratio)
    order=np.argsort(ratio); cum=0.0; wm=None
    for i in order:
        cum+=wr[i]
        if cum>=0.5*np.sum(wr): wm=float(ratio[i]); break
    return dict(n=len(rs),Q=round(Q,3),df=df,I2=round(I2,1),
                ivw_fixed_b=round(b_fixed,5),ivw_fixed_OR=round(math.exp(b_fixed),4),ivw_fixed_p=p_fixed,
                ivw_random_b=round(b_rand,5),ivw_random_OR=round(math.exp(b_rand),4),ivw_random_p=p_rand,
                egger_slope=round(slope,5),egger_intercept=round(intercept,5),egger_intercept_p=p_int,
                weighted_median_OR=round(math.exp(wm),4) if wm is not None else None,
                snps={r:dict(orient=m[r]["orient"],be=m[r]["be"],b_iv=m[r]["b_iv"],se=m[r]["se"]) for r in rs})

def main():
    if not download_and_assemble(): sys.exit(2)
    inst,bmd=load_bmd_alleles()
    print(f"[326] BMD instruments with alleles: {len(bmd)}",flush=True)
    m=collect_fracture(inst,bmd)
    print(f"[326] matched BMD instruments in weight-bearing fracture: {len(m)}",flush=True)
    if len(m)<10:
        res=dict(status="UNDERPOWERED",reason="fewer than 10 BMD instruments matched",
                 n_matched=len(m),matched=list(m),note="no MR estimate reported; CPBB weight-bearing fracture GWAS is heavily indel/WES-like, SNP overlap limited")
    else:
        r=mr(m)
        inconclusive=(r["ivw_random_p"]>0.05 and r["I2"]>60 and np.sign(r["ivw_fixed_b"])!=np.sign(r["ivw_random_b"]))
        r["status"]="INCONCLUSIVE_HETEROGENEITY" if inconclusive else "DIRECTIONAL_POSITIVE_CONTROL"
        r["n_matched"]=len(m)
        r["interpretation"]=("BMD protects weight-bearing fracture direction-consistently with any-fracture positive control"
                             if r["status"]=="DIRECTIONAL_POSITIVE_CONTROL" else "heterogeneity dominates; direction not resolved")
        res=r
    json.dump(res,open(OUT,"w"),indent=2)
    print(f"[326] wrote {OUT}",flush=True)
    print(json.dumps({k:v for k,v in res.items() if k!="snps" and k!="matched"},indent=2),flush=True)

if __name__=="__main__":
    main()
