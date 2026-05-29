import anndata as ad

adata = ad.concat({
    f.name[:-5]: sc.read_h5ad(f)
    for f in tqdm(pathlib.Path("/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/adata/mapped_pathology_data").glob("*.h5ad"))
}, join="outer", label="folder")
adata.write_h5ad("/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/adata/mapped_pathology_data/concatenated_full_1APR2025.h5ad")