def plot_all_chips_grid_expression(
    adata,
    genes=None,
    cell_types=None,
    grid_size=100,
    mode="normal",
    score_mode="product",
    celltype_key="subclass",
    L_gene=None,
    R_gene=None,
    outline_json_dir=None,
    library_id_key="orig.ident",
    cmap="magma_r"
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import numpy as np
    import pandas as pd
    import os
    import json
    from scipy.sparse import issparse

    def apply_transformation_to_coords(coords, transformation_matrix):
        coords_homogeneous = np.hstack([coords, np.ones((coords.shape[0], 1))])
        matrix = np.array(transformation_matrix).reshape(3, 3)
        transformed = coords_homogeneous @ matrix.T
        return transformed[:, :2]

    all_chip_ids = adata.obs[library_id_key].unique()
    n_chips = len(all_chip_ids)

    # ==== 1. 预计算全局网格范围 ====
    global_xmin, global_xmax = np.inf, -np.inf
    global_ymin, global_ymax = np.inf, -np.inf
    for chip_id in all_chip_ids:
        adata_chip = adata[adata.obs[library_id_key] == chip_id]
        coords = adata_chip.obsm["spatial"]
        x_bin = (coords[:, 0] // grid_size).astype(int)
        y_bin = (coords[:, 1] // grid_size).astype(int)
        global_xmin = min(global_xmin, x_bin.min())
        global_xmax = max(global_xmax, x_bin.max())
        global_ymin = min(global_ymin, y_bin.min())
        global_ymax = max(global_ymax, y_bin.max())

    grid_width = global_xmax - global_xmin + 1
    grid_height = global_ymax - global_ymin + 1
    margin =20
    extent = [
        global_xmin - 0.5,
        global_xmax + 0.5,
        global_ymax + 0.5,
        global_ymin - 0.5
    ]

    fig = plt.figure(figsize=(5 * n_chips, 6), dpi=300)

    for idx, chip_id in enumerate(all_chip_ids):
        adata_chip = adata[adata.obs[library_id_key] == chip_id].copy()
        coords = adata_chip.obsm["spatial"]
        x = coords[:, 0]
        y = coords[:, 1]
        x_bin = (x // grid_size).astype(int)
        y_bin = (y // grid_size).astype(int)
        df = pd.DataFrame({"x": x, "y": y, "x_bin": x_bin, "y_bin": y_bin})

        if mode == "normal":
            assert genes is not None and len(genes) > 0
            for gene in genes:
                assert gene in adata.var_names

            X = adata_chip[:, genes].layers["data"]
            if issparse(X): X = X.toarray()
            for i, gene in enumerate(genes):
                df[gene] = X[:, i]

            if cell_types is not None:
                df[celltype_key] = adata_chip.obs[celltype_key].values
                df = df[df[celltype_key].isin(cell_types)]

            gene_means = df.groupby(["x_bin", "y_bin"])[genes].mean().reset_index()

            if score_mode == "product":
                gene_means["score"] = gene_means[genes].prod(axis=1)
            elif score_mode == "sum":
                gene_means["score"] = gene_means[genes].sum(axis=1)
            elif score_mode == "mean":
                gene_means["score"] = gene_means[genes].mean(axis=1)
            else:
                raise ValueError("score_mode 必须是 'product', 'sum' 或 'mean'")

            grouped = gene_means[["x_bin", "y_bin", "score"]]

        elif mode == "CCI":
            assert L_gene is not None and R_gene is not None and len(cell_types) == 2
            source_type, target_type = cell_types

            for gene in L_gene + R_gene:
                assert gene in adata.var_names

            df[celltype_key] = adata_chip.obs[celltype_key].values
            X = adata_chip[:, L_gene + R_gene].layers["data"]
            if issparse(X): X = X.toarray()

            df_L = df[df[celltype_key] == source_type].copy()
            for i, gene in enumerate(L_gene):
                df_L[gene] = X[df[celltype_key] == source_type, i]

            df_R = df[df[celltype_key] == target_type].copy()
            for j, gene in enumerate(R_gene):
                df_R[gene] = X[df[celltype_key] == target_type, len(L_gene) + j]

            # 在每个 grid 上分别计算每个基因的 mean，再相乘
            grouped_L = df_L.groupby(["x_bin", "y_bin"])[L_gene].mean()
            grouped_L["L_expr"] = grouped_L.prod(axis=1)

            grouped_R = df_R.groupby(["x_bin", "y_bin"])[R_gene].mean()
            grouped_R["R_expr"] = grouped_R.prod(axis=1)

            # 合并后相乘
            merged = pd.merge(grouped_L["L_expr"], grouped_R["R_expr"], left_index=True, right_index=True)
            merged["score"] = merged["L_expr"] * merged["R_expr"]
            merged = merged.reset_index()
            grouped = merged[["x_bin", "y_bin", "score"]]
        else:
            raise ValueError("mode 必须是 'normal' 或 'CCI'")
                # ==== 填充统一 score_grid ====
        score_grid = np.full((grid_height, grid_width), np.nan)
        for _, row in grouped.iterrows():
            xi = int(row["x_bin"] - global_xmin)
            yi = int(row["y_bin"] - global_ymin)
            score_grid[yi, xi] = row["score"]

        base_cmap = plt.get_cmap(cmap)
        colors = base_cmap(np.linspace(0, 1, 256))
        #colors[0] = [0.85, 0.85, 0.85, 0.4] # 灰色代表真实为0
        custom_cmap = ListedColormap(colors)
        
        # 只对 score > 0 计算 vmax（避免被 0 拉低）
        positive_scores = score_grid[score_grid > 0]
        vmin = 0
        vmax = np.nanpercentile(positive_scores, 99) if positive_scores.size > 0 else 1
        
        # 指定 Normalize 范围，让 0 映射到 colors[0]
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        # ==== 用 add_axes 精确分布每个子图 ====
        ax = fig.add_axes([0.05 + idx * (0.9 / n_chips), 0.1, 0.9 / n_chips, 0.8])
        im = ax.imshow(
            score_grid,
            cmap=custom_cmap,
            origin="upper",
            interpolation='none',
            extent=extent,
            norm=norm
        )

        ax.set_xlim(global_xmin - 0.5 - margin*1, global_xmax + 0.5 + margin*1)
        ax.set_ylim(global_ymax + 0.5 + margin*5, global_ymin - 0.5 - margin*1.5)
        # ax.set_xlim(global_xmin - 0.5 - margin*0.25, global_xmax + 0.5 + margin*0.25)
        # ax.set_ylim(global_ymax + 0.5 + margin*0.2, global_ymin - 0.5 - margin*0.2)
        ax.set_aspect("equal")
        ax.axis('off')
        ax.set_title(f"{chip_id}", fontsize=12)

        # ==== 轮廓线 ====
        if outline_json_dir is not None:
            json_file = os.path.join(
                f"/sdd/datasets/motor_datasets_107c/mq{outline_json_dir}_cellbins/ROI/outline/MQ{outline_json_dir}-{chip_id}-P0.json"
                #f"/sdd/bgi/wangzilu/00.motorcortex_stduio/motorcortex_proj/41.final_CellChat_3/mq{outline_json_dir}_outline/Mq{outline_json_dir}-{chip_id}-P0.json"
            )
            if os.path.exists(json_file):
                with open(json_file, 'r') as f:
                    data = json.load(f)
                matrix = data['corr_para']['transformation_matrix']
                for line in data['lines']:
                    pts = np.array(line['countours'])
                    pts_trans = apply_transformation_to_coords(pts, matrix)
                    ax.plot(
                        pts_trans[:, 0] / grid_size - 0.5,
                        pts_trans[:, 1] / grid_size - 0.5,
                        color='black',
                        lw=0.5,
                        alpha=0.5
                    )

    plt.show()
    return fig
