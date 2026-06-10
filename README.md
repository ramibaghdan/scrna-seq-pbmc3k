# Single-cell RNA-seq: PBMC 3k (scanpy)

A clustering and cell-type annotation analysis of ~2,700 peripheral blood
mononuclear cells (PBMCs) from a healthy donor, the standard 10x Genomics
benchmark dataset, using [scanpy](https://scanpy.readthedocs.io/).

## What the analysis does

```
load  ->  QC + filter  ->  normalize + log  ->  highly variable genes
      ->  PCA  ->  neighbor graph  ->  UMAP  ->  Leiden clustering
      ->  marker genes  ->  annotate cell types
```

The result is a UMAP where clusters are labeled with immune cell types
(T cells, B cells, NK cells, monocytes, dendritic cells, megakaryocytes),
identified from canonical marker genes.

## Repo layout

| Path | What it is |
|---|---|
| `src/pbmc3k_analysis.py` | The full analysis, written with `# %%` cell markers |
| `environment.yml` / `requirements.txt` | Reproducible environment |
| `figures/` | Saved plots (QC, UMAP, markers, annotated clusters) |
| `notebooks/` | The rendered `.ipynb` (generated from the `.py` via jupytext) |

## Setup + run

```bash
# Option A: conda
conda env create -f environment.yml
conda activate scrna

# Option B: venv + pip
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Run as a script (saves figures to figures/)
python src/pbmc3k_analysis.py

# Or turn it into a notebook with outputs, then open it
jupytext --to notebook --execute src/pbmc3k_analysis.py -o notebooks/pbmc3k_analysis.ipynb
jupyter lab notebooks/pbmc3k_analysis.ipynb
```

## Key markers used for annotation

| Marker | Cell type |
|---|---|
| CD3D | T cells |
| CD8A | Cytotoxic (CD8) T cells |
| MS4A1 | B cells |
| NKG7 | NK cells |
| CD14, LYZ | CD14+ monocytes |
| FCGR3A | CD16+ (FCGR3A) monocytes |
| PPBP | Megakaryocytes / platelets |

## Notes / talking points

- Leiden cluster numbers are not fixed across runs; the cluster-to-cell-type
  mapping in the script must be checked against the marker-gene plots and
  adjusted. (Doing this by hand is the point: it is how you reason about identity.)
- Natural extensions: cell-cycle scoring, batch integration across donors
  (e.g. Harmony/scVI), or differential expression between conditions.

> Built as a self-directed project to work hands-on with the standard single-cell
> RNA-seq analysis workflow.
