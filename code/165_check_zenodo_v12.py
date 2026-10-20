# -*- coding: utf-8 -*-
"""确认 Zenodo v1.2 版本是否已 mint (概念 DOI 10.5281/zenodo.23400860)。

带重试(沙箱代理对 zenodo.org 会间歇性 4xx/断开) + 多源交叉核对。
"""
import json, time, urllib.request, urllib.error

HDR = {"User-Agent": "zenodo-v12-check", "Accept": "application/json"}


def get(url, tries=8, wait=3.0):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url)
            for k, v in HDR.items():
                req.add_header(k, v)
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            time.sleep(wait)
    print("   FAIL:", str(last)[:110])
    return None


def show(h, tag=""):
    md = h.get("metadata", {})
    files = h.get("files") or []
    tot = sum(f.get("size", 0) for f in files) / 1048576
    keys = ", ".join(f.get("key") for f in files)
    print(f"  id={h.get('id')}  doi={h.get('doi')}  version={md.get('version','?')}  created={h.get('created')}")
    print(f"     conceptrecid={h.get('conceptrecid')}  files={len(files)} ({tot:.2f} MB)")
    print(f"     asset: {keys[:100]}")
    return h.get("created"), h.get("doi")


seen = {}

# 源 1: 版本列表 (挂在任意已知版本记录上)
print("=== source 1: /versions on record 23400861 ===")
d = get("https://zenodo.org/api/records/23400861/versions?size=50")
if d:
    hits = d.get("hits", {}).get("hits", [])
    print("   total:", d.get("hits", {}).get("total"), "| returned:", len(hits))
    for h in hits:
        seen[h.get("id")] = h
        show(h)

# 源 2: 精确搜索 (稳健版: 无花哨参数)
print("\n=== source 2: search by repo name ===")
d = get("https://zenodo.org/api/records?q=OA-BMD-fracture-MR&size=50")
if d:
    for h in d.get("hits", {}).get("hits", []):
        seen[h.get("id")] = h
        show(h)

# 源 3: 已知 v1.1 id, 再补一次
print("\n=== source 3: direct fetch of known v1.1 record 23505928 ===")
h = get("https://zenodo.org/api/records/23505928")
if h:
    seen[h.get("id")] = h
    show(h)

print("\n" + "=" * 78)
print("ALL RECORDS (sorted by created):")
print("=" * 78)
for rid in sorted(seen, key=lambda r: str(seen[r].get("created", ""))):
    h = seen[rid]
    md = h.get("metadata", {})
    files = h.get("files") or []
    tot = sum(f.get("size", 0) for f in files) / 1048576
    print(f"  created={h.get('created')}  id={rid}  doi={h.get('doi')}  version={md.get('version','?')}"
          f"  [{len(files)} files, {tot:.2f} MB]")
    print(f"     {(h.get('links') or {}).get('self_html')}")

today_n = len([r for r in seen if str(seen[r].get("created", "")).startswith("2026-10-2")])
print()
print("本日(2026-10-2x)记录数:", today_n)
if today_n >= 2:
    print("VERDICT: v1.2 已 mint ✅ (概念下已有 %d 条本日记录)" % today_n)
elif today_n == 1:
    print("VERDICT: 仅 v1.1 已 mint, v1.2 仍在排队 ⏳")
else:
    print("VERDICT: 尚未出现本日记录 ⏳")
print("NOTE: Zenodo 稿件应引用概念 DOI 10.5281/zenodo.23400860 (恒指向最新版本)")
