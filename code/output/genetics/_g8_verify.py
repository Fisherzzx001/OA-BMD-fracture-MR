# G8 验证：新 13x13 rg_matrix 与 (a) 旧 12x12 逐格一致，(b) fracture 行与 rg_fracture.csv 一致，(c) 与 rg_matrix_13.csv 一致
import pandas as pd, numpy as np

new13 = pd.read_csv('rg_matrix.csv', index_col=0)
old12 = pd.read_csv('_pre_g8_backup/rg_matrix.csv', index_col=0)

blk = new13.loc[old12.index, old12.columns]
d1 = (blk - old12).abs().max().max()
print('(a) max |delta| 12x12 block vs _pre_g8_backup/rg_matrix.csv :', d1)

fr = pd.read_csv('rg_fracture.csv').set_index('trait')['rg']
row = new13.loc['fracture']
mx = 0.0
for t in fr.index:
    d = abs(row[t] - fr[t])
    mx = max(mx, d)
    print(f'    {t:16s} new={row[t]:+.6f} code10={fr[t]:+.6f} d={row[t]-fr[t]:+.2e}')
print('(b) max |delta| fracture row vs rg_fracture.csv :', mx)

m13 = pd.read_csv('_pre_g8_backup/rg_matrix_13.csv', index_col=0)
d3 = (new13.loc[m13.index, m13.columns] - m13).abs().max().max()
print('(c) max |delta| vs _pre_g8_backup/rg_matrix_13.csv :', d3)

# 结构核验
print('shape:', new13.shape, '| diag==1:', bool((np.diag(new13.values) == 1).all()),
      '| symmetric:', bool((new13.values == new13.values.T).all() if new13.values.dtype != object else False))
print('columns:', list(new13.columns))
