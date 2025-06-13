import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.sparse import issparse

def plot_grid_expression_score(
    adata,
    genes=None,
    cell_types=None,
    grid_size=100,
    mode="normal",  # "normal" or "CCI"
    score_mode="product",  # for "normal" mode: "sum", "mean", or "product"
    celltype_key="celltype",
    L_gene=None,
    R_gene=None
):
    assert "spatial" in adata.obsm, "adata.obsm 中必须包含 'spatial' 空间坐标"
    coords = adata.obsm["spatial"]
    x = coords[:, 0]
    y = coords[:, 1]

    x_bin = (x // grid_size).astype(int)
    y_bin = (y // grid_size).astype(int)

    # 所有细胞格子（不含真正空格子）
    all_df = pd.DataFrame({
        "x_bin": x_bin,
        "y_bin": y_bin
    })
    all_counts = all_df.groupby(["x_bin", "y_bin"]).size().reset_index(name="total_cell_count")

    # ✅ 构建全体可能格子（包含真正空格子）
    x_range = np.arange(x_bin.min(), x_bin.max() + 1)
    y_range = np.arange(y_bin.min(), y_bin.max() + 1)
    all_possible_coords = pd.MultiIndex.from_product([x_range, y_range], names=["x_bin", "y_bin"]).to_frame(index=False)

    # 合并所有格子，缺失 total_cell_count 设为 0
    all_counts = pd.merge(all_possible_coords, all_counts, how="left", on=["x_bin", "y_bin"])
    all_counts["total_cell_count"] = all_counts["total_cell_count"].fillna(0)

    df = pd.DataFrame({
        "x": x,
        "y": y,
        "x_bin": x_bin,
        "y_bin": y_bin
    })

    if mode == "normal":
        assert genes is not None and len(genes) > 0, "normal 模式下必须传入 genes"
        for gene in genes:
            assert gene in adata.var_names, f"{gene} 不在 adata.var_names 中"

        X = adata[:, genes].layers["data"]
        if issparse(X):
            X = X.toarray()

        for i, gene in enumerate(genes):
            df[gene] = X[:, i]

        if cell_types is not None:
            assert celltype_key in adata.obs.columns, f"adata.obs 中缺少 '{celltype_key}' 列"
            df[celltype_key] = adata.obs[celltype_key].values
            df = df[df[celltype_key].isin(cell_types)]

        if score_mode == "product":
            df["score"] = df[genes].prod(axis=1)
        elif score_mode == "sum":
            df["score"] = df[genes].sum(axis=1)
        elif score_mode == "mean":
            df["score"] = df[genes].mean(axis=1)
        else:
            raise ValueError("score_mode 必须是 'product', 'sum' 或 'mean'")

        df["cell_count"] = 1
        grouped = df.groupby(["x_bin", "y_bin"]).agg({
            "score": "sum",
            "cell_count": "count"
        }).reset_index()

    elif mode == "CCI":
        assert L_gene is not None and R_gene is not None, "CCI 模式下必须提供 L_gene 和 R_gene"
        assert len(cell_types) == 2, "CCI 模式下必须传入两个 cell_types"

        for gene in L_gene + R_gene:
            assert gene in adata.var_names, f"{gene} 不在 adata.var_names 中"

        celltype_array = adata.obs[celltype_key].values
        df[celltype_key] = celltype_array

        X_L = adata[:, L_gene].layers["data"]
        X_R = adata[:, R_gene].layers["data"]
        if issparse(X_L): X_L = X_L.toarray()
        if issparse(X_R): X_R = X_R.toarray()

        df["L_expr"] = np.nan
        df["R_expr"] = np.nan

        mask_L = celltype_array == cell_types[0]
        mask_R = celltype_array == cell_types[1]

        if len(L_gene) == 1:
            df.loc[mask_L, "L_expr"] = X_L[mask_L, 0]
        else:
            df.loc[mask_L, "L_expr"] = X_L[mask_L].mean(axis=1)

        if len(R_gene) == 1:
            df.loc[mask_R, "R_expr"] = X_R[mask_R, 0]
        else:
            df.loc[mask_R, "R_expr"] = X_R[mask_R].mean(axis=1)

        # 每个格子内两类细胞分别求均值，然后乘积作为得分
        L_group = df[mask_L].groupby(["x_bin", "y_bin"])["L_expr"].mean().reset_index()
        R_group = df[mask_R].groupby(["x_bin", "y_bin"])["R_expr"].mean().reset_index()

        merged_group = pd.merge(L_group, R_group, how="inner", on=["x_bin", "y_bin"])
        merged_group["score"] = merged_group["L_expr"] * merged_group["R_expr"]
        grouped = merged_group[["x_bin", "y_bin", "score"]].copy()
        grouped["cell_count"] = 1  # 占位

    else:
        raise ValueError("mode 必须是 'normal' 或 'CCI'")

    # 合并格子，判断空格子
    merged = pd.merge(all_counts, grouped, how="left", on=["x_bin", "y_bin"])
    merged["score"] = merged["score"].fillna(0)
    merged["cell_count"] = merged["cell_count"].fillna(0)

    max_radius = 50 // grid_size
    is_empty_dict = {
        (row["x_bin"], row["y_bin"]): row["total_cell_count"] == 0
        for _, row in merged.iterrows()
    }

    confirmed_empty_set = set()
    for coord, is_empty in is_empty_dict.items():
        if not is_empty:
            continue
        x0, y0 = coord
        all_neighbors_empty = True
        for dx in range(-max_radius, max_radius + 1):
            for dy in range(-max_radius, max_radius + 1):
                if dx == 0 and dy == 0:
                    continue
                neighbor = (x0 + dx, y0 + dy)
                if neighbor in is_empty_dict and not is_empty_dict[neighbor]:
                    all_neighbors_empty = False
                    break
            if not all_neighbors_empty:
                break
        if all_neighbors_empty:
            confirmed_empty_set.add(coord)

    x_unique = np.arange(x_bin.min(), x_bin.max() + 1)
    y_unique = np.arange(y_bin.min(), y_bin.max() + 1)
    x_map = {x_val: i for i, x_val in enumerate(x_unique)}
    y_map = {y_val: i for i, y_val in enumerate(y_unique)}

    score_grid = np.full((len(y_unique), len(x_unique)), np.nan)

    for _, row in merged.iterrows():
        xi = x_map[row["x_bin"]]
        yi = y_map[row["y_bin"]]
        coord = (row["x_bin"], row["y_bin"])
        if coord in confirmed_empty_set:
            score_grid[yi, xi] = np.nan
        else:
            score_grid[yi, xi] = row["score"]

    # 可视化
    cmap = plt.get_cmap("magma_r")
    colors = cmap(np.linspace(0, 1, 256))
    colors[0] = [0.7, 0.7, 0.7, 1.0]
    custom_cmap = ListedColormap(colors)

    non_zero_scores = score_grid[score_grid > 0]
    if non_zero_scores.size == 0:
        print("⚠️ 所有格子中得分均为 0，跳过颜色映射自动缩放，使用默认 vmin=0, vmax=1")
        vmin, vmax = 0, 1
    else:
        vmin = np.nanmin(non_zero_scores)
        vmax = np.nanmax(score_grid)

    norm = plt.Normalize(vmin=0, vmax=vmax)

    plt.figure(figsize=(10, 8))
    im = plt.imshow(score_grid, cmap=custom_cmap, origin="upper", interpolation='none', norm=norm)
    plt.colorbar(im, label=f"{mode} score")
    plt.title(f"Grid-based spatial {mode} score (grid size = {grid_size})")
    plt.xlabel("Grid X")
    plt.ylabel("Grid Y")
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.show()

    return score_grid
