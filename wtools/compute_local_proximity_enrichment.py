import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm
from collections import defaultdict

def compute_local_proximity_enrichment(
    adata, 
    radius=50, 
    jitter_range=100,
    n_permutations=100,
    seed=42
):
    np.random.seed(seed)
    results = {}

    for area in adata.obs['region'].unique():
        print(f"\nProcessing Area: {area}")
        sub_adata = adata[adata.obs['region'] == area].copy()
        # 去除为 'nan', 'NaN', 'None', 'null' 等无效字符串
        invalid_values = ['nan', 'NaN', 'none', 'None', 'null', 'NULL']
        sub_adata = sub_adata[~sub_adata.obs['celltype'].isin(invalid_values)].copy()
        coords = sub_adata.obs[['tx', 'ty']].values
        labels = sub_adata.obs['celltype'].values.astype(str)
        unique_celltypes = np.unique(labels)
        n_cells = len(coords)

        # Real counts
        tree_real = cKDTree(coords)
        neighbors_real = tree_real.query_ball_point(coords, r=radius)

        real_count = pd.DataFrame(
            0, index=unique_celltypes, columns=unique_celltypes, dtype=int
        )

        for i, nbrs in enumerate(neighbors_real):
            src = labels[i]
            for j in nbrs:
                if i == j: continue
                tgt = labels[j]
                real_count.loc[src, tgt] += 1

        # Null distributions using localized jitter
        null_distributions = defaultdict(lambda: [])

        for p in tqdm(range(n_permutations), desc=f"Permuting {area}"):
            # 局部扰动
            jittered_coords = coords + np.random.uniform(
                low=-jitter_range, high=jitter_range, size=coords.shape
            )
            tree_perm = cKDTree(jittered_coords)
            neighbors_perm = tree_perm.query_ball_point(jittered_coords, r=radius)

            perm_count = pd.DataFrame(
                0, index=unique_celltypes, columns=unique_celltypes, dtype=int
            )

            for i, nbrs in enumerate(neighbors_perm):
                src = labels[i]
                for j in nbrs:
                    if i == j: continue
                    tgt = labels[j]
                    perm_count.loc[src, tgt] += 1

            # 存储每次打乱的计数
            for src in unique_celltypes:
                for tgt in unique_celltypes:
                    null_distributions[(src, tgt)].append(perm_count.loc[src, tgt])

        # 计算p值
        pval_df = pd.DataFrame(index=unique_celltypes, columns=unique_celltypes, dtype=float)

        for src in unique_celltypes:
            for tgt in unique_celltypes:
                real = real_count.loc[src, tgt]
                null = np.array(null_distributions[(src, tgt)])
                pval = (np.sum(null >= real) + 1) / (n_permutations + 1)  # empirical p
                pval_df.loc[src, tgt] = pval

        results[area] = pval_df

    return results
