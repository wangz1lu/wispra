def plot_genes_grid_expression(
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
    cmap="magma_r",
    x_margin_factor_left=5,
    x_margin_factor_right=5,
    y_margin_factor_top=5,
    y_margin_factor_bottom=5
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import numpy as np
    import pandas as pd
    import os
    import json
    import glob
    from scipy.sparse import issparse

    def apply_transformation_to_coords(coords, transformation_matrix):
        coords_homogeneous = np.hstack([coords, np.ones((coords.shape[0], 1))])
        matrix = np.array(transformation_matrix).reshape(3, 3)
        transformed = coords_homogeneous @ matrix.T
        return transformed[:, :2]
    def find_outline_json(outline_json_dir, chip_id):
        pattern = os.path.join(outline_json_dir, f"*{chip_id}*.json")
        matched_files = glob.glob(pattern)
        if not matched_files:
            raise FileNotFoundError(f"No JSON file containing '{chip_id}' found in {outline_json_dir}")
        return matched_files[0]
        
    all_chip_ids = adata.obs[library_id_key].unique()
    n_chips = len(all_chip_ids)


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


            grouped_L = df_L.groupby(["x_bin", "y_bin"])[L_gene].mean()
            grouped_L["L_expr"] = grouped_L.prod(axis=1)

            grouped_R = df_R.groupby(["x_bin", "y_bin"])[R_gene].mean()
            grouped_R["R_expr"] = grouped_R.prod(axis=1)

            merged = pd.merge(grouped_L["L_expr"], grouped_R["R_expr"], left_index=True, right_index=True)
            merged["score"] = merged["L_expr"] * merged["R_expr"]
            merged = merged.reset_index()
            grouped = merged[["x_bin", "y_bin", "score"]]
        else:
            raise ValueError("mode 必须是 'normal' 或 'CCI'")

        score_grid = np.full((grid_height, grid_width), np.nan)
        for _, row in grouped.iterrows():
            xi = int(row["x_bin"] - global_xmin)
            yi = int(row["y_bin"] - global_ymin)
            score_grid[yi, xi] = row["score"]

        base_cmap = plt.get_cmap(cmap)
        colors = base_cmap(np.linspace(0, 1, 256))

        custom_cmap = ListedColormap(colors)
        

        positive_scores = score_grid[score_grid > 0]
        vmin = 0
        vmax = np.nanpercentile(positive_scores, 99) if positive_scores.size > 0 else 1
        
 
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

        ax = fig.add_axes([0.05 + idx * (0.9 / n_chips), 0.1, 0.9 / n_chips, 0.8])
        im = ax.imshow(
            score_grid,
            cmap=custom_cmap,
            origin="upper",
            interpolation='none',
            extent=extent,
            norm=norm
        )
        ax.set_xlim(global_xmin - 0.5 - margin * x_margin_factor_left,
            global_xmax + 0.5 + margin * x_margin_factor_right)
        ax.set_ylim(global_ymax + 0.5 + margin * y_margin_factor_top,
            global_ymin - 0.5 - margin * y_margin_factor_bottom)

        ax.set_aspect("equal")
        ax.axis('off')
        ax.set_title(f"{chip_id}", fontsize=12)

        # ==== 轮廓线 ====
        if outline_json_dir is not None:
            json_file = find_outline_json(outline_json_dir, chip_id)
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

def plot_colocalization_celltypes(
    adata,
    cell_types,
    celltype_key="subclass",
    library_id_key="library_id",
    point_size=1,
    figsize=(5, 5),
    celltype_colors=["red", "blue"],
    background_color="lightgray",
    distance_threshold=400,
    outline_json_dir=None,
    grid_size=100,
    x_margin_factor_left=500,
    x_margin_factor_right=500,
    y_margin_factor_top=500,
    y_margin_factor_bottom=500,
    show_legend=True,
    invert_y=False
):

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.spatial.distance import cdist
    import os
    import json
    import glob
    def apply_transformation_to_coords(coords, transformation_matrix):
        coords_homogeneous = np.hstack([coords, np.ones((coords.shape[0], 1))])
        matrix = np.array(transformation_matrix).reshape(3, 3)
        transformed = coords_homogeneous @ matrix.T
        return transformed[:, :2]    
    def find_outline_json(outline_json_dir, chip_id):
        pattern = os.path.join(outline_json_dir, f"*{chip_id}*.json")
        matched_files = glob.glob(pattern)
        if not matched_files:
            raise FileNotFoundError(f"No JSON file containing '{chip_id}' found in {outline_json_dir}")
        return matched_files[0]
    chips = adata.obs[library_id_key].unique()
    n_chips = len(chips)


    global_xmin, global_xmax = np.inf, -np.inf
    global_ymin, global_ymax = np.inf, -np.inf

    for chip_id in chips:
        coords = adata[adata.obs[library_id_key] == chip_id].obsm["spatial"]
        global_xmin = min(global_xmin, coords[:, 0].min())
        global_xmax = max(global_xmax, coords[:, 0].max())
        global_ymin = min(global_ymin, coords[:, 1].min())
        global_ymax = max(global_ymax, coords[:, 1].max())

    margin = 80  
    extent = [
        global_xmin - margin,
        global_xmax + margin,
        global_ymax + margin,
        global_ymin - margin
    ]

    fig = plt.figure(figsize=(figsize[0] * n_chips, figsize[1]), dpi=300)

    for i, chip_id in enumerate(chips):
        ax = fig.add_axes([0.05 + i * (0.9 / n_chips), 0.1, 0.9 / n_chips, 0.8])

        adata_chip = adata[adata.obs[library_id_key] == chip_id]
        coords = adata_chip.obsm["spatial"]
        ct = adata_chip.obs[celltype_key]

        idx1 = ct == cell_types[0]
        idx2 = ct == cell_types[1]
        coords1 = coords[idx1]
        coords2 = coords[idx2]
            
        dist_matrix = cdist(coords1, coords2)
        close_pairs = np.where(dist_matrix < distance_threshold)
        idx1_keep = np.unique(np.where(idx1)[0][close_pairs[0]])
        idx2_keep = np.unique(np.where(idx2)[0][close_pairs[1]])

        keep_mask = np.zeros(len(adata_chip), dtype=bool)
        keep_mask[idx1_keep] = True
        keep_mask[idx2_keep] = True

        other_mask = ~(keep_mask & (idx1 | idx2))
        

        coords_plot = coords.copy()
        if invert_y:
            coords_plot[:, 1] = -coords_plot[:, 1]
        

        ax.scatter(
            coords_plot[other_mask, 0], coords_plot[other_mask, 1],
            c=background_color, s=point_size * 1, alpha=1, edgecolors='none'
        )

        ax.scatter(
            coords_plot[idx1_keep, 0], coords_plot[idx1_keep, 1],
            c=celltype_colors[0], s=point_size * 10, label=cell_types[0], alpha=1, edgecolors='none'
        )
        ax.scatter(
            coords_plot[idx2_keep, 0], coords_plot[idx2_keep, 1],
            c=celltype_colors[1], s=point_size * 10, label=cell_types[1], alpha=1, edgecolors='none'
        )



        json_file = find_outline_json(outline_json_dir, chip_id)
        if os.path.exists(json_file):
            with open(json_file, 'r') as f:
                data = json.load(f)
            matrix = data['corr_para']['transformation_matrix']
            for line in data['lines']:
                pts = np.array(line['countours'])
                pts_trans = apply_transformation_to_coords(pts, matrix)
                if invert_y:
                    pts_trans[:, 1] = -pts_trans[:, 1]
                ax.plot(
                    pts_trans[:, 0],
                    pts_trans[:, 1],
                    color='black',
                    lw=0.5,
                    alpha=0.5
                )
        ax.set_xlim(global_xmin - 5 - margin*x_margin_factor_left, global_xmax + 5 + margin*x_margin_factor_right)
        ax.set_ylim(global_ymax + 5 + margin*y_margin_factor_top, global_ymin - 5 - margin*y_margin_factor_bottom)

        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.axis("off")
        ax.set_title(f"{chip_id}", fontsize=10)


    if show_legend:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=celltype_colors[0], edgecolor='none', label=cell_types[0]),
            Patch(facecolor=celltype_colors[1], edgecolor='none', label=cell_types[1]),
            Patch(facecolor=background_color, edgecolor='none', label="Other")
        ]

        fig.subplots_adjust(bottom=0.15)  
        fig.legend(
            handles=legend_elements,
            loc='lower center',
            bbox_to_anchor=(0.5, 0.2),  
            ncol=3,
            frameon=False,
            handlelength=1.2,
            fontsize=10
        )

    plt.show()
    return fig
