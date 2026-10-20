# -*- coding: utf-8 -*-
"""
G8 追加审计 5（终局）：拦截 code/10 生产运行实际传给 compute_gencov 的参数
=======================================================================
审计 4 显示：用完全复刻的输入（桥接 ld 字典 + 原始文件桥接 Zdf + code/10 参数顺序）
调用 compute_gencov，得到 gc=0.0029546012 / rg=+0.118958 —— 与 code/10 生产运行的
+0.148117 不符。所有输入量（ld 字典、Zdf、ell、Ns、h2、截距）均已证明一致。

本脚本猴子补丁 m08.compute_gencov / m08.compute_h2，记录 code/10.main() 实际
传入的参数与 compute_h2 的返回值，首次 gencov 调用后即中止，随后与审计 4 的
复刻值逐项对比，定位分歧的精确位置。
"""
import os
import sys
import importlib
import numpy as np

CODE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE)

m10 = importlib.import_module("10_fracture_row")
m08 = importlib.import_module("08_genetic_correlation")

h2_calls = []
gc_calls = []

orig_h2 = m08.compute_h2
orig_gc = m08.compute_gencov


def logged_h2(z, ell, N, M, n_iter=3, n_blocks=200):
    r = orig_h2(z, ell, N, M, n_iter, n_blocks)
    h2_calls.append(dict(N=N, M=M, n_iter=n_iter, n_blocks=n_blocks,
                         z_head=z[:5].copy(), z_tail=z[-5:].copy(),
                         ell_head=ell[:5].copy(),
                         h2=r[0], intercept=r[2]))
    return r


def logged_gc(z1, z2, ell, N1, N2, M, h1=None, h2=None, i_h1=1.0, i_h2=1.0,
              n_iter=3, n_blocks=200):
    gc_calls.append(dict(
        z1_head=z1[:5].copy(), z2_head=z2[:5].copy(),
        z1_tail=z1[-5:].copy(), z2_tail=z2[-5:].copy(),
        ell_head=ell[:5].copy(), ell_tail=ell[-5:].copy(),
        N1=N1, N2=N2, M=M, h1=h1, h2=h2, i_h1=i_h1, i_h2=i_h2,
        n_iter=n_iter, n_blocks=n_blocks,
        len_y=len(z1)))
    raise SystemExit("captured first gencov call")


m08.compute_h2 = logged_h2
m08.compute_gencov = logged_gc

try:
    m10.main()
except SystemExit:
    pass

print(f"\n=== compute_h2 调用 {len(h2_calls)} 次 ===")
for i, c in enumerate(h2_calls):
    print(f"  [{i}] N={c['N']:,.1f} M={c['M']:,.0f} blocks={c['n_blocks']} "
          f"h2={c['h2']:.12f} intercept={c['intercept']:.9f}")
    print(f"       z_head={np.array2string(c['z_head'], precision=4)}")
    print(f"       ell_head={np.array2string(c['ell_head'], precision=2)}")

print(f"\n=== compute_gencov 调用 {len(gc_calls)} 次 ===")
for c in gc_calls:
    print(f"  N1={c['N1']:,.1f} N2={c['N2']:,.1f} M={c['M']:,.0f} "
          f"h1={c['h1']:.12f} h2={c['h2']:.12f} i_h1={c['i_h1']:.9f} i_h2={c['i_h2']:.9f} "
          f"n_iter={c['n_iter']} n_blocks={c['n_blocks']} len_y={c['len_y']:,}")
    print(f"  z1_head={np.array2string(c['z1_head'], precision=4)}")
    print(f"  z2_head={np.array2string(c['z2_head'], precision=4)}")
    print(f"  z1_tail={np.array2string(c['z1_tail'], precision=4)}")
    print(f"  z2_tail={np.array2string(c['z2_tail'], precision=4)}")
    print(f"  ell_head={np.array2string(c['ell_head'], precision=2)}")
    print(f"  ell_tail={np.array2string(c['ell_tail'], precision=2)}")

np.save(os.path.join(CODE, "output", "genetics", "_g8_gc_z1.npy"), gc_calls[0]["z1_head"] if gc_calls else np.array([]))
