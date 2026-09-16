"""
07_download_ldsc_ref.py — 从 Zenodo(7768714) 分块下载 S-LDSC 参考文件
=========================================================================
只用 urllib（绕过 Windows schannel 的证书吊销问题），分块写入 + 断点续传 + 重试。
目标：
  - 1000G_Phase3_ldscores.tgz      (ref LD 计分，若含 EUR QC l2 则用)
  - 1000G_Phase3_weights_hm3_no_MHC.tgz (回归权重 w-ld)
  - 1000G_Phase3_frq.tgz           (MAF，munge 用)
  - 1000G_Phase3_baselineLD_v2.2_ldscores.tgz (若上面不含 EUR QC l2，则 baselineLD 内含)
"""
import os
import sys
import time
import tarfile
import urllib.request

REC = "7768714"
BASE = f"https://zenodo.org/api/records/{REC}/files"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref")
os.makedirs(OUT, exist_ok=True)

WANT = [
    "1000G_Phase3_ldscores.tgz",
    "1000G_Phase3_weights_hm3_no_MHC.tgz",
    "1000G_Phase3_frq.tgz",
    "1000G_Phase3_baselineLD_v2.2_ldscores.tgz",
]
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}


def download(key, max_retry=5):
    url = f"{BASE}/{key}/content"
    dest = os.path.join(OUT, key)
    for attempt in range(max_retry):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                total = int(r.headers.get("Content-Length", 0))
                done = 0
                with open(dest, "wb") as f:
                    while True:
                        buf = r.read(1 << 20)
                        if not buf:
                            break
                        f.write(buf)
                        done += len(buf)
                print(f"  ✓ {key}  {done/1e6:.1f} MB")
                return True
        except Exception as e:
            print(f"  ✗ {key} attempt {attempt+1}: {type(e).__name__} {e}")
            if os.path.exists(dest):
                os.remove(dest)
            time.sleep(3)
    return False


def list_tgz(key):
    dest = os.path.join(OUT, key)
    if not os.path.exists(dest):
        return
    try:
        with tarfile.open(dest, "r:gz") as t:
            names = t.getnames()
        eur = [n for n in names if "EUR.QC" in n and n.endswith(".l2.ldscore.gz")]
        w = [n for n in names if "weights" in n and n.endswith(".l2.ldscore.gz")]
        m = [n for n in names if n.endswith(".l2.M")]
        print(f"\n[{key}] 共 {len(names)} 项 | EUR.QC l2={len(eur)} | weights l2={len(w)} | .l2.M={len(m)}")
        for n in (eur[:3] + w[:2] + m[:2]):
            print("    ", n)
    except Exception as e:
        print(f"  列举失败 {key}: {e}")


if __name__ == "__main__":
    print("▶ 下载 S-LDSC 参考文件")
    for k in WANT:
        download(k)
    print("\n▶ 检查内容")
    for k in WANT:
        list_tgz(k)
    print("\n完成。若 baselineLD 内含 EUR.QC l2，则用它作 --ref-ld-chr。")
