import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .pl import plot_genes_grid_expression

def compute_interaction_score(
    adata,
    interaction_dict,
    output_dir,
    plot_name_prefix="result",
    grid_size=400,
    score_mode="product",
    cell_types=None,
    celltype_key="SubClass",
    do_scaling=True,
    verbose=True
):
    """
    Compute spatial interaction scores across region/layer combinations and save plots and score matrix.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics AnnData object.
    interaction_dict : dict
        Format: {"interaction_name": [gene1, gene2], ...}
    output_dir : str
        Directory to save results.
    plot_name_prefix : str
        Prefix for output file names (PDF/CSV).
    grid_size : int
        Size of spatial grid.
    score_mode : str
        Scoring mode, e.g., "product".
    cell_types : list or None
        Subset of cell types to include.
    celltype_key : str
        Key in `adata.obs` for cell type annotation.
    do_scaling : bool
        Whether to apply z-score scaling (not used here).
    verbose : bool
        Whether to print progress.

    Returns
    -------
    heatmap_wide : pd.DataFrame
        Interaction score matrix (wide format).
    skipped_pairs : list
        List of interaction names that were skipped due to error or missing data.
    """
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = os.path.join(output_dir, f"{plot_name_prefix}_LRpairs.pdf")
    csv_path = os.path.join(output_dir, f"{plot_name_prefix}_score_matrix.csv")

    heatmap_matrix_long = []
    skipped_pairs = []

    with PdfPages(pdf_path) as pdf:
        for name, gene_list in interaction_dict.items():
            if verbose:
                print(f"Plotting {name}: {gene_list}")
            try:
                score_grid, fig, df_grid = plot_genes_grid_expression(
                    adata=adata,
                    grid_size=grid_size,
                    mode="normal",
                    genes=gene_list,
                    score_mode=score_mode,
                    cell_types=cell_types,
                    celltype_key=celltype_key
                )

                pdf.savefig(fig)
                plt.close(fig)

                non_zero = df_grid["score"][df_grid["score"] > 0]
                score_max = np.nanpercentile(non_zero, 99.5) if len(non_zero) > 0 else 1.0
                df_grid["score_capped"] = df_grid["score"].clip(upper=score_max)

                if "region" in df_grid.columns and "layer" in df_grid.columns:
                    df_grouped = (
                        df_grid.groupby(["region", "layer"])["score_capped"]
                        .mean()
                        .reset_index()
                    )
                    df_grouped["interaction"] = name
                    heatmap_matrix_long.append(df_grouped)
                else:
                    if verbose:
                        print(f"{name} missing region/layer columns, skipping")

            except Exception as e:
                if verbose:
                    print(f"Skipping {name} due to error: {e}")
                skipped_pairs.append(name)

    if heatmap_matrix_long:
        heatmap_df = pd.concat(heatmap_matrix_long, axis=0)
        heatmap_df["region_layer"] = heatmap_df["region"] + "_" + heatmap_df["layer"]
        heatmap_wide = heatmap_df.pivot(index="interaction", columns="region_layer", values="score_capped")
        heatmap_wide = heatmap_wide.dropna(axis=1)
        heatmap_wide.to_csv(csv_path)
        if verbose:
            print(f"Score matrix saved to: {csv_path}")
    else:
        heatmap_wide = pd.DataFrame()
        if verbose:
            print("No valid interaction scores generated.")

    if verbose:
        print("Plotting complete.")
        print(f"Skipped interactions: {skipped_pairs if skipped_pairs else 'None'}")

    return heatmap_wide, skipped_pairs

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from tqdm import tqdm

def compute_cci_score_sum(
    adata,
    source_type,
    target_type,
    L_gene,
    R_gene,
    celltype_key="subclass",
    spatial_key="spatial",
    library_key="chip",
    grid_size=400
):
    """
    Compute total colocalization score between source and target cell types
    using ligand-receptor product expression on spatial grid.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial and expression data.
    source_type : str
        Cell type label for ligand-expressing source cells.
    target_type : str
        Cell type label for receptor-expressing target cells.
    L_gene : list of str
        List of ligand gene names.
    R_gene : list of str
        List of receptor gene names.
    celltype_key : str, default "subclass"
        Column in `adata.obs` containing cell type labels.
    spatial_key : str, default "spatial"
        Key in `adata.obsm` for spatial coordinates.
    library_key : str, default "chip"
        Column in `adata.obs` identifying different tissue slices or chips.
    grid_size : int, default 400
        Grid bin size (in microns) used to group cells spatially.

    Returns
    -------
    total_score : float
        Total colocalization interaction score across all chips.
    """
    total_score = 0
    for chip_id in adata.obs[library_key].unique():
        adata_chip = adata[adata.obs[library_key] == chip_id].copy()
        coords = adata_chip.obsm[spatial_key]
        x = coords[:, 0]
        y = coords[:, 1]
        x_bin = (x // grid_size).astype(int)
        y_bin = (y // grid_size).astype(int)

        df = pd.DataFrame({"x_bin": x_bin, "y_bin": y_bin})
        df[celltype_key] = adata_chip.obs[celltype_key].values

        genes = L_gene + R_gene
        for gene in genes:
            assert gene in adata.var_names

        X = adata_chip[:, genes].layers["data"]
        if issparse(X):
            X = X.toarray()

        df_L = df[df[celltype_key] == source_type].copy()
        for i, gene in enumerate(L_gene):
            df_L[gene] = X[df[celltype_key] == source_type, i]

        df_R = df[df[celltype_key] == target_type].copy()
        for j, gene in enumerate(R_gene):
            df_R[gene] = X[df[celltype_key] == target_type, len(L_gene) + j]

        if df_L.empty or df_R.empty:
            continue

        grouped_L = df_L.groupby(["x_bin", "y_bin"])[L_gene].mean()
        grouped_L["L_expr"] = grouped_L.prod(axis=1)

        grouped_R = df_R.groupby(["x_bin", "y_bin"])[R_gene].mean()
        grouped_R["R_expr"] = grouped_R.prod(axis=1)

        merged = pd.merge(grouped_L["L_expr"], grouped_R["R_expr"], left_index=True, right_index=True)
        merged["score"] = merged["L_expr"] * merged["R_expr"]

        total_score += merged["score"].sum()

    return total_score


def generate_cci_heatmap(
    adata,
    source_subclasses=None,
    target_subclasses=None,
    subclass_key="subclass",
    celltype_key=None,  
    L_gene=["SLC1A3", "GLS"],
    R_gene=["GRM7"],
    grid_size=400,
    spatial_key="spatial",
    library_key="chip",
    cmap="Reds"
):
    """
    Generate a heatmap of total colocalization interaction scores
    between source and target cell types or subclasses.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial and expression data.
    source_subclasses : list of str
        List of source subclasses (or cell types) to evaluate as ligand producers.
    target_subclasses : list of str
        List of target subclasses (or cell types) to evaluate as receptor producers.
    subclass_key : str, default "subclass"
        Column in `adata.obs` with subclass annotations.
    celltype_key : str or None, default None
        Column in `adata.obs` with finer cell type annotations.
        If None, subclass_key will be used as default.
    L_gene : list of str
        List of ligand gene names.
    R_gene : list of str
        List of receptor gene names.
    grid_size : int, default 400
        Grid bin size (in microns) used to group cells spatially.
    spatial_key : str, default "spatial"
        Key in `adata.obsm` for spatial coordinates.
    library_key : str, default "chip"
        Column in `adata.obs` identifying different tissue slices or chips.
    cmap : str, default "Reds"
        Colormap used for heatmap display.

    Returns
    -------
    heatmap_df : pd.DataFrame
        Heatmap matrix of colocalization interaction scores between each source-target pair.
    """
    if celltype_key is not None:

        df = adata.obs[[subclass_key, celltype_key]]

        source_annots = df[df[subclass_key].isin(source_subclasses)][celltype_key].unique()
        target_annots = df[df[subclass_key].isin(target_subclasses)][celltype_key].unique()

        heatmap_df = pd.DataFrame(index=target_annots, columns=source_annots, dtype=float)

        for t_ann in tqdm(target_annots, desc="Target"):
            for s_ann in source_annots:
                score = compute_cci_score_sum(
                    adata=adata,
                    source_type=s_ann,
                    target_type=t_ann,
                    L_gene=L_gene,
                    R_gene=R_gene,
                    celltype_key=celltype_key,
                    spatial_key=spatial_key,
                    library_key=library_key,
                    grid_size=grid_size
                )
                heatmap_df.at[t_ann, s_ann] = score
    else:
        subclasses = adata.obs[subclass_key]
        source_set = sorted(set(source_subclasses))
        target_set = sorted(set(target_subclasses))

        heatmap_df = pd.DataFrame(index=target_set, columns=source_set, dtype=float)

        for t in tqdm(target_set, desc="Target"):
            for s in source_set:
                score = compute_cci_score_sum(
                    adata=adata,
                    source_type=s,
                    target_type=t,
                    L_gene=L_gene,
                    R_gene=R_gene,
                    celltype_key=subclass_key,  
                    spatial_key=spatial_key,
                    library_key=library_key,
                    grid_size=grid_size
                )
                heatmap_df.at[t, s] = score

    heatmap_df = heatmap_df.fillna(0)

    plt.figure(figsize=(6, 4))
    sns.heatmap(heatmap_df, cmap=cmap, annot=False, fmt=".1f")
    plt.title("CCI Score Sum by Source/Target")
    plt.xlabel("Source Cell Types")
    plt.ylabel("Target Cell Types")
    plt.tight_layout()
    plt.show()

    return heatmap_df
