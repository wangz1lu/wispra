# wispra 空间转录组网格化基因邻近表达可视化工具包

一个旨在优化空间组数据处理和细胞通讯下游分析的包，通过网格化计算策略，优化现有细胞通讯空间可视化绘图体验。

## 功能 (更新中)

* Ligand-Receptor Pairs 空间邻近表达可视化工具
* 任意数量的基因空间邻近表达可视化工具
* 邻近细胞类型空间可视化工具
* 计算、统计以及可视化大量 Ligand-Receptor Pairs 空间分布特征
* 计算、统计以及可视化大量 Source -> Target 细胞对对 Ligand-Receptor Pairs 空间分布特征的贡献
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

## 绘图模块（pl模块）
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
    celltype_key="celltype", # 设置输入细胞类型对应的标签
    grid_size=400, # 设置正方形网格的半径（单位：像素）
    mode="CCI", # CCI 模式
    library_id_key="chip", # 支持多张片子（必选）
    outline_json_dir="/path/to/folder/", # 支持轮廓线（可选）
)

# 邻近细胞类型空间可视化
wp.pl.plot_colocalization_celltypess(
    adata=adata,
    cell_types=["celltype1", "celltype2"],
    celltype_key="celltype", # 设置输入细胞类型对应的标签
    library_id_key="chip", # 支持多张片子
    point_size=0.1, # 点的大小
    celltype_colors=["#CC0000", "#0000CC"], # 设置细胞的颜色
    background_color="#DCDCDC", # 设置背景点的颜色
    distance_threshold=200, # 设置两个细胞类型之间的距离
    outline_json_dir="/path/to/folder/", # 支持轮廓线（可选）
    invert_y=True # 反转y轴
)

## 工具模块（tl模块）
# 计算、统计以及可视化大量Ligand-Receptor Pairs空间分布特征
wp.tl.compute_interaction_score(
    adata=adata,
    interaction_dict, # LR pairs字典，键是pathway名称，值是对应的LR pairs
    output_dir,
    plot_name_prefix="result", #结果文件命名
    grid_size=400, # 设置正方形网格的半径（单位：像素）
)

# 计算、统计以及可视化大量 Source -> Target 细胞对对 Ligand-Receptor Pairs 空间分布特征的贡献
wp.tl.generate_cci_heatmap(
    adata=adata,
    source_subclasses=["subclass1","subclass2","subclass3"], # 指定source细胞类型
    target_subclasses=["subclass4","subclass5","subclass6"], # 指定target细胞类型
    subclass_key="subclass", # 设置输入细胞类型对应的标签
    celltype_key=None, # 如果提供精度更高的细胞类型注释，会自动计算精度高的细胞对之间的分数（可选）
    L_gene=["gene1"], # 支持配体复合物["gene1","gene2"]
    R_gene=["gene2"], # 支持受体复合物["gene3","gene4"]
    grid_size=400, # 设置正方形网格的半径（单位：像素）
    library_key="chip", # 支持多张片子（必选）
)
```
注：所有pl函数均支持参数调整绘图显示范围，这在绘制轮廓线时很有帮助
    x_margin_factor_left=500
    x_margin_factor_right=500
    y_margin_factor_top=500
    y_margin_factor_bottom=500
    
