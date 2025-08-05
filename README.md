# Wispra 空间组处理和数据可视化工具包

一个旨在优化空间组数据处理和细胞通讯下游分析的包，优化现有细胞通讯空间可视化绘图体验。

## 功能 (更新中)

* 任意数量的基因空间邻近表达可视化工具
* Ligand-Receptor Pairs 空间邻近表达可视化工具
* 邻近细胞类型空间可视化工具
* 计算、统计以及可视化大量Ligand-Receptor Pairs空间分布特征
* 计算、统计以及可视化大量Source -> Target细胞对对Ligand-Receptor Pairs空间分布特征的贡献
* spatial数据处理独家工具包（待更新）
  
## 安装

```bash
pip install wispra
```
## 更新
```bash
pip install --upgrade wispra
```
## 用法
```python
import wispra as wp

# Ligand-Receptor Pairs 空间邻近表达可视化
wp.pl.plot_genes_grid_expression(
    adata=adata,
    genes=["gene1","gene2"],
    grid_size=400,
    mode="normal",
    library_id_key="chip" #多张片子
)

# Cell-to-Cell interaction 事件可视化
wp.pl.plot_genes_grid_expression(
    adata=adata,
    L_genes=["gene1"], #支持配体复合物
    R_genes=["gene2"], #支持受体复合物
    cell_types=["celltype1","celltype2"]，
    celltype_key="celltype",
    grid_size=400,
    mode="CCI",
    library_id_key="chip" #多张片子
)

# 邻近细胞类型空间可视化工具

```
