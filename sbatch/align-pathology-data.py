import sys
import pathlib
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.spatial
from tqdm import tqdm
from matplotlib.path import Path

folder_key = sys.argv[1]

example = list(pathlib.Path("/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/raw/xenium/").glob(f"*{folder_key}*"))[0]

lipid_droplets = pd.read_csv(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/locations/lipid-droplet/{folder_key}.csv").drop("Unnamed: 0", axis=1).mul(0.2125)
plin2 = pd.read_csv(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/locations/plin2/{folder_key}.csv").drop("Unnamed: 0", axis=1).mul(0.2125)
oil_red_o = pd.read_csv(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/locations/oil-red-o/{folder_key}.csv").drop("Unnamed: 0", axis=1).mul(0.2125)
amyloid = pd.read_csv(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/locations/amyloid/{folder_key}.csv").drop("Unnamed: 0", axis=1)

def read_xenium_data(xenium_path: pathlib.Path):
    adata = sc.read_10x_h5(xenium_path / "cell_feature_matrix.h5")
    adata.obs = pd.read_csv(xenium_path / "cells.csv.gz").set_index("cell_id")
    adata.uns["cell_boundary"] = {cell_id: cell_boundaries[["vertex_x", "vertex_y"]].to_numpy().reshape(-1) for cell_id, cell_boundaries in tqdm(pd.read_csv(xenium_path / "cell_boundaries.csv.gz").groupby("cell_id"))}

    adata.layers["transcript"] = pd.read_csv(xenium_path / "transcripts.csv.gz") \
        .eval("near_nucleus = nucleus_distance <= 5") \
        .pivot_table(
            index="cell_id",
            columns="feature_name",
            values="near_nucleus",
            aggfunc="sum",
            fill_value=0
        ) \
        .reindex(index=adata.obs_names, columns=adata.var_names, fill_value=0).values

    return adata

def assign_distance_measurement(adata: ad.AnnData, measurement: pd.DataFrame, measurement_name: str):
    # Setup cKD Tree
    ckd_tree = scipy.spatial.cKDTree(measurement[["x_centroid", "y_centroid"]])

    # Query measurement 
    distance, index = ckd_tree.query(adata.obs[["x_centroid", "y_centroid"]])
    adata.obs[measurement_name] = distance
    return adata

def assign_stain_measurement(adata: ad.AnnData, measurement: pd.DataFrame, measurement_name: str):
    # Setup cKD Tree
    ckd_tree = scipy.spatial.cKDTree(adata.obs[["x_centroid", "y_centroid"]])

    # Query measurement 
    distance, index = ckd_tree.query(measurement[["x_centroid", "y_centroid"]])
    measurement["minimal_distance"] = distance
    measurement["cell_id"] = [adata.obs.index[ix] for ix in index]

    adata.obs[measurement_name] = 0
    adata.obs[measurement_name.replace("area", "count")] = 0
    for row in tqdm(measurement.itertuples(), total = measurement.shape[0]):
        if Path(adata.uns["cell_boundary"][row.cell_id].reshape(-1, 2)).contains_point((row.x_centroid, row.y_centroid)):
            adata.obs.loc[row.cell_id, measurement_name] += row.area
            adata.obs.loc[row.cell_id, measurement_name.replace("area", "count")] += 1
    
    return adata

adata = read_xenium_data(example)

adata = assign_stain_measurement(adata, lipid_droplets, "lipid_droplet_area")
adata = assign_stain_measurement(adata, plin2, "plin2_area")
adata = assign_stain_measurement(adata, oil_red_o, "oil_red_o_area")
adata = assign_distance_measurement(adata, amyloid, "distance_to_nearest_amyloid")

adata.write_h5ad(f"/oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/data/processed/adata/mapped_pathology_data/{folder_key}.h5ad")