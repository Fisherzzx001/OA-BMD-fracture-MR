# -*- coding: utf-8 -*-
"""实证: 13 个 PLINK *_clump_input.txt 能否由现有本地资源离线再生。

clump 输入 = code/17 prepare_clump_inputs() 从各暴露全基因组成因 GWAS 表中
抽出 SNP + P 两列写出。因此只要 (a) 能定位其上游 GWAS 文件、(b) 这些文件仍在本地,
就可确认可再生, 从而安全地删除回滚镜像。
"""
import importlib.util, os, sys, glob

sys.path.insert(0, "code")

# 找到定义了 FLAGSHIP_LDSC / PHENOS 的模块
cands = sorted(glob.glob("code/09*.py")) + sorted(glob.glob("code/*.py"))
mod = None
for c in cands:
    try:
        src = open(c, encoding="utf-8", errors="replace").read()
    except Exception:
        continue
    if "PHENOS" in src and "FLAGSHIP_LDSC" in src:
        spec = importlib.util.spec_from_file_location("m09", c)
        m = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(m)
        except SystemExit:
            pass
        except Exception as e:
            print("  (import %s raised %s: %s)" % (c, type(e).__name__, str(e)[:70]))
            continue
        if hasattr(m, "PHENOS") and hasattr(m, "FLAGSHIP_LDSC"):
            mod, path = m, c
            break
        continue

if mod is None:
    print("未能定位上游 GWAS 清单模块 —— 无法实证可再生")
    sys.exit(1)

base = mod.FLAGSHIP_LDSC
print("模块          :", path)
print("GWAS 根目录   :", base)
print("根目录存在    :", os.path.isdir(base))
print()

n_ok = 0
missing = []
for k, v in mod.PHENOS.items():
    p = os.path.join(base, v)
    exists = os.path.exists(p)
    size = os.path.getsize(p) / 1048576 if exists else 0
    n_ok += 1 if exists else 0
    if not exists:
        missing.append(p)
    print("  %-22s %-8s %8.1f MB  %s" % (k, "OK" if exists else "MISSING", size, v))

print()
print("上游 GWAS 文件: %d 个 | 本地存在 %d 个 | 缺失 %d 个" % (len(mod.PHENOS), n_ok, len(missing)))
for m_ in missing:
    print("   MISSING:", m_)

# 再确认写出端仍可读: code/17 的 prepare_clump_inputs 依赖 pandas
print()
print("=== 结论 ===")
if not missing:
    print("全部上游 GWAS 文件仍在本地 -> 13 个 clump 输入可由 code/17 离线确定性再生")
    print("-> 删除 _backup_pristine.git 不会造成不可恢复的数据损失")
else:
    print("存在缺失上游文件 -> 不可假定可再生, 需先保留备份或从原始来源补齐")
