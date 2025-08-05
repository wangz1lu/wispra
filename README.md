# Wispra 空间转录组处理和数据可视化工具包

一个旨在优化空间组数据处理和细胞通讯下游分析的包，优化现有细胞通讯空间可视化绘图体验。

## 功能 (更新中)

* Ligand-Receptor Pairs 空间邻近表达可视化工具
* 任意数量的基因空间邻近表达可视化工具
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
    genes=["gene1","gene2"], # 支持单/多基因，即任意数量的基因空间邻近表达可视化
    grid_size=400, #设置正方形网格的半径（单位：像素）
    mode="normal", # 普通模式
    library_id_key="chip",  # 支持多张片子（必选）
    outline_json_dir="/path/to/folder/", # 支持轮廓线（可选）
) 

# Cell-to-Cell interaction 事件可视化
wp.pl.plot_genes_grid_expression(
    adata=adata,
    L_genes=["gene1"], # 支持配体复合物["gene1","gene2"]
    R_genes=["gene2"], # 支持受体复合物["gene3","gene4"]
    cell_types=["celltype1","celltype2"]，
    celltype_key="celltype",
    grid_size=400, #设置正方形网格的半径（单位：像素）
    mode="CCI", # CCI 模式
    library_id_key="chip", # 支持多张片子（必选）
    outline_json_dir="/path/to/folder/", # 支持轮廓线（可选）
)

# 邻近细胞类型空间可视化工具
wp.pl.plot_colocalization_celltypess(
    adata=adata,
    cell_types=["celltype1", "celltype2"],
    celltype_key="celltype",
    library_id_key="chip", # 支持多张片子
    point_size=0.1, # 点的大小
    celltype_colors=["#CC0000", "#0000CC"], # 设置细胞的颜色
    background_color="#DCDCDC", # 设置背景点的颜色
    distance_threshold=200,
    outline_json_dir="/path/to/folder/", # 支持轮廓线（可选）
    x_margin_factor_left=500,
    x_margin_factor_right=500,
    y_margin_factor_top=-400,
    y_margin_factor_bottom=1500,
    invert_y=True # 反转y轴
)
```
