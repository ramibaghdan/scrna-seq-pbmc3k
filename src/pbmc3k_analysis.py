# ---------------------------------------------------------------------------
# Single-cell RNA-seq analysis of the PBMC 3k dataset using scanpy.
#
# This file uses `# %%` cell markers, so it runs as a plain script
# (`python src/pbmc3k_analysis.py`) AND opens as notebook cells in VS Code /
# Jupyter. We convert it to a real .ipynb with figures using jupytext.
#
# The workflow is the standard single-cell pipeline:
#   load -> QC + filter -> normalize -> highly variable genes -> PCA
#        -> neighbors -> UMAP -> Leiden clustering -> marker genes -> annotate
# ---------------------------------------------------------------------------

# %% [markdown]
# # PBMC 3k single-cell RNA-seq (scanpy)
# ~2,700 peripheral blood mononuclear cells from a healthy donor (10x Genomics).
# Goal: cluster the cells and identify the immune cell types present.

# %% Imports and settings
import os
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 2  # show pipeline progress
sc.settings.set_figure_params(dpi=80, facecolor="white")

FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

# %% Load data (downloads ~5 MB the first time, then cached locally)
adata = sc.datasets.pbmc3k()
print(adata)  # AnnData: cells (obs) x genes (var)

# %% [markdown]
# ## Quality control
# Real droplets contain one healthy cell. We flag empty/dying cells and doublets
# using: number of genes per cell, total counts, and the fraction of
# mitochondrial reads (high mito % usually means a stressed/dying cell).

# %% Compute QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # mitochondrial genes
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4, multi_panel=True, show=False,
)
plt.savefig(f"{FIGDIR}/01_qc_violin.png", bbox_inches="tight")
plt.close()

# %% Filter low-quality cells and rarely-detected genes
sc.pp.filter_cells(adata, min_genes=200)   # drop cells with too few genes
sc.pp.filter_genes(adata, min_cells=3)     # drop genes seen in <3 cells
adata = adata[adata.obs.n_genes_by_counts < 2500].copy()  # drop likely doublets
adata = adata[adata.obs.pct_counts_mt < 5].copy()         # drop dying cells
print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

# %% [markdown]
# ## Normalization
# Normalize each cell to the same total counts, then log-transform so highly
# expressed genes do not dominate. We keep a copy of the full gene set in `.raw`
# for plotting marker genes later.

# %% Normalize + log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %% Highly variable genes (the informative ones for structure)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(f"{FIGDIR}/02_highly_variable_genes.png", bbox_inches="tight")
plt.close()

adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()

# Regress out technical effects, then scale each gene to unit variance
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# %% [markdown]
# ## Dimensionality reduction + clustering
# PCA compresses the data, a neighbor graph connects similar cells, UMAP lays it
# out in 2D for visualization, and Leiden finds communities (clusters) in the graph.

# %% PCA
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig(f"{FIGDIR}/03_pca_variance.png", bbox_inches="tight")
plt.close()

# %% Neighbor graph + UMAP + Leiden clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)

sc.pl.umap(adata, color=["leiden"], legend_loc="on data", show=False)
plt.savefig(f"{FIGDIR}/04_umap_leiden.png", bbox_inches="tight")
plt.close()

# %% [markdown]
# ## Marker genes
# For each cluster, rank the genes that are most differentially expressed. These
# define the cluster's identity.

# %% Rank marker genes per cluster
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
plt.savefig(f"{FIGDIR}/05_marker_genes.png", bbox_inches="tight")
plt.close()

# %% Canonical PBMC markers projected on the UMAP
# CD3D=T cells, CD8A=cytotoxic T, MS4A1=B cells, NKG7=NK, CD14/LYZ=monocytes,
# FCGR3A=CD16 monocytes, PPBP=megakaryocytes/platelets.
canonical = ["CD3D", "CD8A", "MS4A1", "NKG7", "CD14", "LYZ", "FCGR3A", "PPBP"]
sc.pl.umap(adata, color=canonical, show=False)
plt.savefig(f"{FIGDIR}/06_umap_canonical_markers.png", bbox_inches="tight")
plt.close()

# %% [markdown]
# ## Annotate clusters
# Map each Leiden cluster to a cell type using the markers above.
# IMPORTANT: cluster numbers depend on the run. Inspect figures 05/06 and adjust
# this mapping so each label matches the markers its cluster expresses.

# %% Label clusters (EDIT after inspecting the marker plots)
cluster_to_celltype = {
    "0": "CD4 T",
    "1": "CD14+ Monocytes",
    "2": "B",
    "3": "CD8 T",
    "4": "NK",
    "5": "FCGR3A+ Monocytes",
    "6": "Dendritic",
    "7": "Megakaryocytes",
}
adata.obs["cell_type"] = (
    adata.obs["leiden"].map(cluster_to_celltype).astype("category")
)
sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False)
plt.savefig(f"{FIGDIR}/07_umap_celltypes.png", bbox_inches="tight")
plt.close()

# %% Save the processed object for reuse
adata.write("pbmc3k_processed.h5ad")
print("Done. Figures in figures/, processed data in pbmc3k_processed.h5ad")
